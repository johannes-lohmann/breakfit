# breakfit

This code accompanies the paper "Objective extraction and analysis of statistical features of Dansgaard–Oeschger events" by Johannes Lohmann and Peter D. Ditlevsen.

The purpose is to fit a continuous piecewise-linear model to a time series with repeated saw-tooth shaped cycles, such as the Dansgaard-Oeschger events of the last glacial period.
This is achieved by an interative procedure, which is described in more detail in the paper. 

To use the code, compile the cython file (.pyx) with the command

python setup.py build_ext --inplace

and then execute the main .py file.

When using the code, please cite the original publication:

Johannes Lohmann and Peter D. Ditlevsen, "Objective extraction and analysis of statistical features of Dansgaard–Oeschger events",
Clim. Past, 15, 1771–1792, 2019,
https://doi.org/10.5194/cp-15-1771-2019
