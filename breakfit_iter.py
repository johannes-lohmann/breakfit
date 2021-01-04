import numpy as np
import matplotlib.pyplot as pl
import time as tm
import breakfit_cython_ext as fit_breakmodel
from scipy.optimize import minimize, basinhopping

def MAIN():

        # initial guesses for the DO warming onsets, and the DO stadial onsets.
        DO_events = np.loadtxt('DO_times_RAS14_coarse.txt',usecols=(1,),unpack=True)
        DO_up = DO_events[9:-3:2]
        DO_down = DO_events[10:-2:2]
        DO_down = np.insert(DO_down,1,85800.)
        
        # record to perform the piecewise-linear fit
        time, record = np.loadtxt('greenland_d18o_stack_11k.txt')

        slopes_u =[];slopes_d =[];up_dur =[];down_dur =[]
        iterations = 45

        start_time = tm.time()
        for j in range(iterations):

                up_idx = [find_nearest(time, x) for x in DO_up]
                down_idx = [find_nearest(time, x) for x in DO_down]

                stad_levels = [np.mean(record[up_idx[i]:down_idx[i]][::-1]) for i in range(1,len(up_idx))]
                # individual DO cycles according to the break points of the current iteration
                transitions = [record[up_idx[i+1]:down_idx[i]][::-1]-stad_levels[i-1] for i in range(1,len(up_idx)-1)]
                # time points for the individual segments
                timings = [[time[up_idx[i]:down_idx[i]],time[down_idx[i+1]:up_idx[i]],time[up_idx[i+1]:down_idx[i+1]]] for i in range(1,len(up_idx)-1)]
                # durations of individual segments
                durations = [[len(timings[i-1][0]),len(timings[i-1][1]),len(timings[i-1][2])] for i in range(1,len(up_idx)-1)]

                results = []
                for i in range(1,len(up_idx)-1):
                        result, value = fit_optimize(transitions,stad_levels,timings,slopes_u,slopes_d,up_dur,down_dur,i,j)
                        results.append(result)

                np.save('stack_iter%s'%j,results)
                
                print(results)

                #Update breakpoints etc as input to next iteration.
                interstadials = [results [i][5]-results [i][0] for i in range(len(results))]

                # first stadial duration
                stadial0 = results[0][0]
                # remaining stadial durations
                stadials = [durations[i][0] + (results[i][0] - durations[i][0]) + (durations[i-1][0]+ durations[i-1][1] - results[i-1][5]) for i in range(1,len(up_idx)-2)]
                stadials = np.append(stadial0,stadials)
                slopes_u = [results[i][1] for i in range(len(results))]
                slopes_d = [results[i][3] for i in range(len(results))]
                up_dur = [results [i][2]-results [i][0] for i in range(len(results))]
                down_dur = [results [i][5]-results [i][4] for i in range(len(results))]

                DO_up_new = [DO_down[1] - (np.sum(stadials[0:i+1]) + np.sum(interstadials[0:i])) for i in range(len(stadials))][::1]
                DO_down_new = [DO_down[1] - (np.sum(stadials[0:i]) + np.sum(interstadials[0:i])) for i in range(len(stadials)+1)][::1]
                DO_up_new = np.append(DO_up[0], DO_up_new)
                DO_up_new = np.append(DO_up_new,DO_up[-1])
                DO_down_new = np.append(DO_down[0], DO_down_new)
                DO_down_new = np.append(DO_down_new,DO_down[-1])

                DO_down = DO_down_new
                DO_up = DO_up_new

        print('simulation time (min) :', (tm.time()-start_time)/60.)

def fit_optimize(transitions,stad_levels,timings,slopes_u,slopes_d,up_dur,down_dur,choose,it):

        offset = stad_levels[choose-1]
        constant = stad_levels[choose]
        transition = transitions[choose-1]
        length = len(transition)
        print(length)
        durations = [len(timings[choose-1][0]),len(timings[choose-1][1]),len(timings[choose-1][2])]
        print(durations)
        if it>0:
                print(up_dur[choose-1],down_dur[choose-1])

        def f(p):
                # some constraints hard-coded into the goal function:
                if p[2]-2<=p[0]:#constraint: up-break before top-break: LIMIT to 2 years
                        return 1.9 +(p[0]-p[2]+2)*0.1
                elif p[4]<=p[2]:#constraint: top-break before down-break
                        return 1.9+(p[2]-p[4])*0.1
                elif p[5]-2<=p[4]:#constraint: down-break before stadial-break: LIMIT to 2 years
                        return 1.9+(p[4]-p[5]+2)*0.1
                elif np.isnan(p[0]):
                        return 10.
                else:
                        return fit_breakmodel.calc_error(transition,np.append(p,constant-offset))


        if it==0:
                # guesses for the parameters at initial iteration
                p0 = np.asarray([durations[0]-10,0.056*1000,durations[0]+50,-0.00146*10000,durations[0]+durations[1]-30,durations[0]+durations[1]+10])
        else:
                p0 = np.asarray([durations[0],slopes_u[choose-1],durations[0]+up_dur[choose-1],slopes_d[choose-1],durations[0]+durations[1]-down_dur[choose-1],durations[0]+durations[1]])

        cons = []
        ###CONSTRAINTS for the constrained optimization.
        b1_b2 = {'type': 'ineq', 'fun' : lambda x: np.array([x[2]-x[0]-2])}
        b2_b3 = {'type': 'ineq', 'fun' : lambda x: np.array([x[4]-x[2]])}
        b3_b4 = {'type': 'ineq', 'fun' : lambda x: np.array([x[5]-x[4]-2])}
        no_below = {'type': 'ineq', 'fun' : lambda x: np.array([((x[2]-x[0])*x[1]/1000.+(x[4]-x[2])*x[3]/10000.)-(constant-offset)])}
        gradual_longer = {'type': 'ineq', 'fun' : lambda x: np.array([x[4]-x[2]-2.*(x[5]-x[4])])}
        drop_steeper = {'type': 'ineq', 'fun' : lambda x: np.array([2*x[3]/10000.*(x[5]-x[4])+x[1]/1000.*(x[2]-x[0])+x[3]/10000.*(x[4]-x[2])-(constant-offset)])}

        b1_low = {'type': 'ineq', 'fun' : lambda x: np.array([x[0]-20])}
        b2_low = {'type': 'ineq', 'fun' : lambda x: np.array([x[2]-20])}
        b4_high = {'type': 'ineq', 'fun' : lambda x: np.array([length-20-x[4]])}
        b5_high = {'type': 'ineq', 'fun' : lambda x: np.array([length-20-x[5]])}
        s1_low = {'type': 'ineq', 'fun' : lambda x: np.array([x[1]-0.02*1000])}
        s1_high = {'type': 'ineq', 'fun' : lambda x: np.array([1.5*1000-x[1]])}
        s2_low = {'type': 'ineq', 'fun' : lambda x: np.array([x[3]+0.3*10000])}
        s2_high = {'type': 'ineq', 'fun' : lambda x: np.array([-0.0001*10000-x[3]])}

        cons.extend((b1_b2,b2_b3,b3_b4,no_below,gradual_longer,drop_steeper,b1_low,b2_low,b4_high,b5_high,s1_low,s1_high,s2_low,s2_high))
        

        class MyTakeStep(object):
                def __init__(self, stepsize=0.5):
                        self.stepsize = stepsize
                def __call__(self, x):
                        # random steps of the basinhopping proposal distribution
                        k0 = x[0] + np.random.normal(loc=0.0,scale=15.)
                        k1 = x[1] + np.random.normal(loc=0.0,scale=0.004*1000)
                        k2 = x[2] + np.random.normal(loc=0.0,scale=15.)
                        k3 = x[3] + np.random.normal(loc=0.0,scale=0.0015*10000)
                        k4 = x[4] + np.random.normal(loc=0.0,scale=15.)
                        k5 = x[5] + np.random.normal(loc=0.0,scale=15.)
                        count1 = 0; count2 = 0
                        # make sure that basinhopping proposal fulfills constraints
                        while ((k2-2<k0) or (k1/1000.<0.02) or (k1/1000.>1.5) or (k0<20) or (k2<20) or ((k2-k0)*k1/1000.<1.8) or ((k2-k0)*k1/1000.>7.5)) and (count1<100000):
                                k0 = x[0] + np.random.normal(loc=0.0,scale=15.)
                                k1 = x[1] + np.random.normal(loc=0.0,scale=0.004*1000)
                                k2 = x[2] + np.random.normal(loc=0.0,scale=15.)
                                count1+=1
                        # make sure that basinhopping proposal fulfills constraints
                        while ((k4<k2) or (k5-2<k4) or (k3/10000.<-0.3) or (k3/10000.>-0.0001) or (k4<durations[0]+0.5*durations[1]) or (k5<20) or (k5>length-10) or ((k1/1000.*(k2-k0)+k3/10000.*(k4-k2))<(constant-offset)) or ((k4-k2)<2.*(k5-k4)) or (((k2-k0)*k1/1000.+(k4-k2)*k3/10000.-(constant-offset))<-2.*k3/10000.*(k5-k4))) and (count2<100000):
                                k3 = x[3] + np.random.normal(loc=0.0,scale=0.0015*10000)
                                k4 = x[4] + np.random.normal(loc=0.0,scale=15.)
                                k5 = x[5] + np.random.normal(loc=0.0,scale=15.)
                                count2+=1
                        x = np.asarray([k0,k1,k2,k3,k4,k5])
                        return x

        mytakestep = MyTakeStep()

        def print_fun(x, f, accepted):
                print("%s" % (f))

        minimizer_kwargs = {"method":"COBYLA", "jac":False, "constraints":cons, "tol":0.0001, "options":{"maxiter":1000, "rhobeg": 0.1} }#L-BFGS-B, ,SLSQP; "bounds":bnds,
        res = basinhopping(f, p0, minimizer_kwargs=minimizer_kwargs, T=0.01,stepsize=5.0, niter=5000, take_step=mytakestep)#, callback=print_fun; niter=5000
        print(res)

        return res.x, res.fun

###find index of timeseries where value is closest to given.
def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx


if __name__ == '__main__':
        MAIN()
