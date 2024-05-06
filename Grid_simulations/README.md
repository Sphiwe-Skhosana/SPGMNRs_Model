This folder contains all the results obtained from evaluating the sensitivity of the proposed MB-EM algorithm on the set of grid points.
We considered the following choices of the number and/or set of grid points:
* $N=50$ chosen uniformly over the domain of $x$
* $N=n$ and set $\mathcal{U}=x$
* $N=100$ chosen uniformly over the domain of $x$
* $N=500$ chosen uniformly over the domain of $x$
  
We generated $200$ samples of size $n=1000$ from the model 

$f(Y|X=x)=0.65\mathcal{N}\\{Y|1-\text{cos}(2\pi x),0.09\\}+0.35\mathcal{N}\\{Y|\text{exp}(2x),0.16\\}$

The results are virtually the same (see Performance_n=1000.pdf).
 
