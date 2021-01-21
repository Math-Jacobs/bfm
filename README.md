

# The back-and-forth method in optimal transport

This repository contains the source code used in the paper [A fast approach to optimal transport: The back-and-forth method](https://arxiv.org/pdf/1905.12154.pdf) [1]. The original code was written in C and we provide here a MATLAB wrapper to the C code.



# Documentation
Available here: <https://back-and-forth.netlify.app>.



# Installation

Requirements: FFTW ([download here](http://www.fftw.org/)), MATLAB.

Download the C MEX file `w2.c` [here](lalala) or clone this GitHub repository.

Compilation: in a MATLAB session run
```matlab
mex -O CFLAGS="\$CFLAGS -std=c99" -lfftw3 -lm w2.c 
```
This will produce a MEX function `w2` that you can use in MATLAB. You may need to use flags `-I` and `-L` to link to the FFTW3 library, e.g. `mex -O CFLAGS="\$CFLAGS -std=c++11" w2.c -lfftw3 -I/usr/local/include`. See [this page](https://www.mathworks.com/help/matlab/matlab_external/build-an-executable-mex-file.html) for more information on how to compile MEX files. 



# Usage

In a MATLAB session, run the command
```matlab
[phi, psi] = w2(mu, nu, numIters, sigma);
```

Input:

* `mu` and `nu` are two arrays of nonnegative values which sum up to the same value.
* `numIters` is the total number of iterations.
* `sigma` is the initial step size of the gradient ascent iterations.

Output:

* `phi` and `psi` are arrays corresponding to the Kantorovich potentials. 





# References


[1] Matt Jacobs and Flavien Léger. [A fast approach to optimal transport: The back-and-forth method](https://arxiv.org/pdf/1905.12154.pdf). *Numerische Mathematik* (2020): 1-32.




