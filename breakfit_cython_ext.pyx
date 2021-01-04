import numpy as np 
cimport numpy as np
cimport cython
ctypedef np.float64_t DT

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.binding(True)
cpdef calc_error(np.ndarray[DT, ndim=1, negative_indices=False, mode='c'] data, np.ndarray[DT, ndim=1, negative_indices=False, mode='c'] params):
    cdef int i
    cdef double errors=0
    cdef double result
    cdef int N = len(data)
    cdef np.ndarray[DT, ndim=1, negative_indices=False, mode='c'] model,model1,model2,model3,model4,model5

    model1=np.zeros(int(params[0]))
    model2=np.linspace(0.,(params[2]-params[0])*params[1]/1000.,(params[2]-params[0]))
    model3=np.linspace((params[2]-params[0])*params[1]/1000.,(params[2]-params[0])*params[1]/1000.+(params[4]-params[2])*params[3]/10000.,(params[4]-params[2]))
    model4=np.linspace((params[2]-params[0])*params[1]/1000.+(params[4]-params[2])*params[3]/10000.,params[6],(params[5]-params[4]))
    model5=np.asarray(int(N-params[5])*[params[6]])
    model = np.concatenate((model1,model2,model3,model4,model5))

    for i in xrange(N):
        errors += (data[i]-model[i])**2
    result = np.sqrt(errors/N)
    return result


