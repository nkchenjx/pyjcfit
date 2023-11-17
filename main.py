# Coded by Jixin Chen @ Ohio University, Department of Chemistry and Biochemistry
# First coded in MATLAB on 2016/04/05
# Converted to python 3.12 2023/11/07. Developed and tested on software versions Python 3.12, Windows 11, Anaconda 3 11/2023.

# MIT License
# Copyright (c) 2023 Jixin Chen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# Anyone uses or modifies pyjcfit function please keep the above copyright information.

import numpy as np
import pandas as pd
from matplotlib import pyplot
from scipy.optimize import curve_fit
from pyjcfit import pyjcfit
import time

# ----------- An example objective function and main function  ---------------:
# def objective(x, para):
#     return para[0] * x + para[1]
def objective(x, para):
    return np.add(np.multiply(para[0],x), para[1])


if __name__ == '__main__':

    # # load input variables
    # data = pd.read_csv("data.csv")
    # # print(data)
    # x_value = data.x
    # y_value = data.y
    # #   print(y_value)

    # or simulate a set of data
    n = 1000
    x_value = np.arange(0, n, 1)
    para_true = [1.2345, 2,3456]
    y_value = list(np.add(objective(x_value, para_true), np.random.normal(0, 10, len(x_value))))


    # initial guess of parameters and bounds
    para_guess = [1, 1]  # initial guessed parameters for the objective function
    bounds = {'ub': [1000, 1000], 'lb': [-1000, -1000]}  # upper bounds (ub) and lower bounds (lb) ordered the same as parameters
    # check if bounds and para_guess are well-ordered.
    d1 = np.subtract(bounds['ub'], para_guess)
    d2 = np.subtract(para_guess, bounds['lb'])
    d3 = np.multiply(d1, d2)
    for i in range(len(d3)):
        if d3[i] < 0:
            print(" para_guess[%s] and its bounds out of order." % i)

    # set searching options
    option = {'maxiteration': 100, 'precision': 0.00001, 'exp_step': 0.5, 'convgtest': 1E-100, 'show_error_surface': True}
    # maxiteration is the maximum searching iteration.
    # precision defines the significant figures. It is the smallest numerical search step of each paramter. e.g. paraguess of previous iteration = 10 and precision = 0.01, then searching step is 0.1 for this iteration and this parameter, i.e. precision = 0.01 is 2 sig fig.
    # exp_step, searching step size +-para*precision*(2^exp_step)^n where n is 1, 2, 3,...
    # convgtest is the minimum difference of sum(residual**2) between searching steps to judge convergence.

    #  --------- fitting starts
    start_time = time.time()
    fit_results = pyjcfit(objective, x_value, y_value, para_guess, bounds, option)
    print('\n---- elapsed time for pyjcfit = %f seconds ----' % (time.time() - start_time))
    para = fit_results['para']
    para_hist = fit_results['para_hist']
    error_hist = fit_results['error_hist']
    residual = np.subtract(y_value, objective(x_value, para))
    goodness_of_fit = fit_results['gof']

    # ------- fitting ends-----


    # ------ plot results
    print('\nfitted para = ', para)
    print('goodness of fitting: ', goodness_of_fit)
    print(fit_results['para_hist'])


    x_new = np.arange(0, np.max(x_value), 0.1)
    y_new = objective(x_new, para)
    pyplot.figure()
    pyplot.subplot(1, 2, 1)
    pyplot.scatter(x_value, y_value)
    pyplot.plot(x_new, y_new, '-', color='red')
    pyplot.plot(x_value, residual, '--', color='orange')
    pyplot.title(str('pyjcfit: y =%.5f * x + (%.5f)' % (para[0], para[1])))
    # pyplot.show()

    pyplot.subplot(1, 2, 2)
    pyplot.plot(error_hist, '-', color='red')
    pyplot.title('error history')
    # pyplot.show()

# end of example

    # ------ control function curve_fit running time ----
    print('\n---- comparing function curve_fit from scipy using Levenberg-Marquardt algorithm ----')
    def objectivec(x, a, b):  # for scipy curve_fit
        return a * x + b

    start_time = time.time()
    para, _ = curve_fit(objectivec, x_value, y_value, method='lm') # lm: the Levenberg-Marquardt algorithm through least square residual search.
    print('fitted para', para)
    print('---- elapsed time for curve_fit = %f seconds  ----' % (time.time() - start_time))
    x_new = np.arange(0, np.max(x_value), 0.1)
    y_new = objectivec(x_new, para[0], para[1])
    pyplot.figure()
    pyplot.scatter(x_value, y_value)
    pyplot.plot(x_new, y_new, '-', color='red')
    pyplot.title(str('LM: y =%.5f * x + (%.5f)' % (para[0], para[1])))
    pyplot.show()
