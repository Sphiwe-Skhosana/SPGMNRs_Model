This folder contains all the `R` code used to compare the results when using different interpolation methods. 
* Polynomial interpolation, more specifically Lagrange interpolation, may result in negative values of the variance function. Moreover, due to its numerical instability, it may overshoot and result in large values of the variance function.
* We compared the results based on cubic spline interpolation and linear interpolation and found no difference in the results. The results can be seen in the file named "Performance_n=2000.pdf" for a sample of size $n=2000$. 
