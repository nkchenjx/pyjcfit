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
import time
import pandas as pd
from matplotlib import pyplot
from scipy.optimize import curve_fit
from pyjcfit import pyjcfit
from sklearn.linear_model import LinearRegression

# ----------- An example objective function and main function  ---------------:
# def objective(x, para):
#     return para[0] * x + para[1]


if __name__ == '__main__':

    # # load input variables
    # data = pd.read_csv("data.csv")
    # # print(data)
    # x_value = data.x
    # y_value = data.y
    # #   print(y_value)

    def objective(x, para):
        y = np.add(para[0], np.multiply(para[1], np.exp(np.multiply((-para[2]), x))))
        return y

    # or simulate a set of data
    n = 100
    para_true = (0.2, 1, 10/n)
    x_value = np.arange(0, n, 1)
    # y_value = list(objective(x_value, para_true))
    x_valuen = x_value + np.random.normal(0, 0.00001, len(x_value))
    y_value = list(np.add(objective(x_valuen, para_true), np.random.normal(0, 0.001, len(x_value)))) #SNR about sum of paraTrue(A)/noise



# ------ control function curve_fit running time ----
    print('\n---- control function curve_fit from scipy using Levenberg-Marquardt algorithm ----')

    def exp_decay(x, b, a, k): # for scipy curve_fit
        return b + a * np.exp(x * (-k))


    start_time = time.time()
    para, _ = curve_fit(exp_decay, x_value, y_value, method='lm')
    print('single exponential b a k:', para)
    finish_time = time.time()
    print('fitted para', para)
    print('---- elapsed time for curve_fit = %.9f seconds  ----' % (finish_time - start_time))
    x_new = np.arange(0, np.max(x_value), 0.1)
    y_new = exp_decay(x_new, para[0], para[1], para[2])
    pyplot.figure(100)
    # pyplot.subplot(2, 1, 2)
    pyplot.scatter(x_value, y_value)
    pyplot.plot(x_new, y_new, '-', color='red')
    pyplot.plot(x_value, np.subtract(y_value, exp_decay(x_value, para[0], para[1], para[2])), color='orange')
    pyplot.title(('curve_fit',  para))
    # pyplot.show()





# ------JCFIT starts here initial guess of parameters and bounds
    print('\n---- JCFit starts ---- ')
    para_guess = [1, 1, 1] # initial guessed parameters for the objective function
    bounds = {'ub': (10, 10, 10),
                       'lb': (0, 0, 0)}  # upper bounds (ub) and lower bounds (lb) ordered the same as parameters
    # check if bounds and para_guess are well-ordered.
    d1 = np.subtract(bounds['ub'], para_guess)
    d2 = np.subtract(para_guess, bounds['lb'])
    d3 = np.multiply(d1, d2)
    for i in range(len(d3)):
        if d3[i] < 0:
            print(" para_guess[%s] and its bounds out of order." % i)

    # set searching options
    option = {'maxiteration': 1000, 'precision': 1E-10, 'exp_step': 0.5, 'convgtest': 0, 'score_func': 'L2', 'show_error_surface': False}
    # maxiteration is the maximum searching iteration.
    # precision defines the significant figures. It is the smallest numerical search step of each paramter. e.g. paraguess of previous iteration = 10 and precision = 0.01, then searching step is 0.1 for this iteration and this parameter, i.e. precision = 0.01 is 2 sig fig.
    # exp_step, searching step size +-para*precision*(2^exp_step)^n where n is 1, 2, 3,...
    # convgtest is the minimum difference of sum(residual**2) between searching steps to judge convergence.

    #  --------- fitting starts
    start_time = time.time()
    fit_results = pyjcfit(objective, x_value, y_value, para_guess, bounds, option)
    finish_time = time.time()
    print('\n---- elapsed time for pyjcfit = %.9f seconds ----' % (finish_time - start_time))
    para = fit_results['para']
    para_hist = fit_results['para_hist']
    error_hist = fit_results['error_hist']
    residual = np.subtract(y_value, objective(x_value, para))
    goodness_of_fit = fit_results['gof']

    # ------- fitting ends-----

    # ------ plot results
    print('\nfitted para = ', end = '')
    for i in para:
        print('%.5f' %i, end = ' ')
    print('\nTrue para   = ', end='')
    for i in para_true:
        print('%.5f' % i, end=' ')
    print('\ngoodness of fitting: ', goodness_of_fit)
    # print(fit_results['para_hist'])

    x_new = np.arange(np.min(x_value), np.max(x_value), 0.1)
    y_new = objective(x_new, para)
    pyplot.figure(200)
    pyplot.subplot(2, 1, 1)
    pyplot.scatter(x_value, y_value)
    pyplot.plot(x_new, y_new, '-', color='red')
    pyplot.plot(x_value, residual, '-', color='orange')
    pyplot.title(('pyjcfit: ',  para))
    # pyplot.show()

    pyplot.subplot(2, 1, 2)
    pyplot.plot(error_hist, '-', color='red')
    pyplot.title('error history')


    pyplot.figure()
    pyplot.scatter(x_value, np.log(np.subtract(y_value,para[0])))

    pyplot.show()



# end of example

