# pyjcfit
# Coded by Jixin Chen @ Ohio University, Department of Chemistry and Biochemistry
# The algorithm is ~10x slower than Levenberg–Marquardt searching algorithm, and is suitable for super COMPLICATED model, global fitting, or neural network optimization.
First coded in MATLAB 2014b on 2016/04/05
Converted to python 3.12 on 2023/11/07

Curve fitting for 1D vectors or ND matrix with any given model equations using a random search algorithm. 2D and above data will need to revise the pyjcfit function to allow data loading and treating. Or simply vectorize multidimential data and feed to the function will work without change.

Developed and tested on software versions
Python 3.12
Windows 11
Anaconda 3 11/2023

Basic idea: 
For a given raw data and a math model, use the inital guess of the parameters in the model to generate the guessed data and get the residual between the two sets of data.
The sum of the absolute values of all residuals (L1) or the sum of the square of the residuals (L2) is a function to the parameters for which the LM use an algorithm modified from Gauss-Newton algorithm (https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm) in finding the minimum. While fitnguess will guess a series of possible values within the boundaries of each parameter and take the minimum as the initial guess for the next iteration. The searching space is exponentially distributed along the initial guess to the boundaries e.g. -10 to 10 with initial guess 0 and searching precision ~0.1: [-10, -5, -2, -1, -0.5, -0.2, -0.1, [0], 0.1, 0.2, 0.5, 1, 2, 5, 10]. It does not follow the assumption in LM that these is only one minimum. It is 5-20 times slower than LM in this example, but has a chance to find the "global" minimum in this parameter dimension, which unfortunitly may be different from the global minimum of all parameters. The searching goes through all parameters sequentially in each iteration making it a P searching algorithm (https://en.wikipedia.org/wiki/P_versus_NP_problem). 
This method is modified from the Random Search method especially the Friedman-Savage procedure (https://en.wikipedia.org/wiki/Random_search). (Friedman, M.; Savage, L.J. (1947). Planning experiments seeking maxima, chapter 13 of Techniques of Statistical Analysis, edited by Eisenhart, Hastay, and Wallis. McGraw-Hill Book Co., New York.)

Change your fitting model in the function objective if you use custom functions for a complicated model or global fitting.

cite:
https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.6b05697
Jixin Chen, Joseph R Pyle, Kurt Waldo Sy Piecco, Anatoly B Kolomeisky, Christy F Landes, A Two-Step Method for smFRET Data Analysis, J. Phys. Chem. B, 2016, 120 (29), pp 7128–7132

A semi-exhaustive and brutal searching algorithm using the least-square method to find the best fit of a curve. The least-square method can be changed easily to other residual treatment methods.

The fitting time of each iteration is scaled with the number of parameters of the fitting, the search speed is not very sensitive to the accuracy because of the exponential searching spacing design (not set as an option now and has to be changed manually inside the function). Each parameter is searched one by one within its boundaries. Thus a P^N question is reduced to PN with a cost of losing space coverage.

The code can be changed to fit multiple curves globally by introducing different models for different curves using the same parameters.

An example can be found at:
Juvinch R. Vicente, Ali Rafiei Miandashti, Kurt Sy Piecco, Joseph R. Pyle, Martin E. Kordesch, Jixin Chen*, Single-Particle Organolead Halide Perovskite Photoluminescence as a Probe for Surface Reaction Kinetics. ACS Applied Matierals & Interfaces, 2019, 11(19), 18034-18043.
