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
    def gaussian_peak(x, para): # para: a, b, c in a*exp(-(x-b)^2/2/c^2)
        return para[0]*np.exp(np.divide(np.square(np.subtract(x, para[1])), -2*para[2]**2))


    def data_simulation(x, para):
        y = np.add(np.add(np.multiply(para[0][0], np.square(np.subtract(x, para[0][3]))), np.multiply(para[0][1], np.subtract(x, para[0][3]))), para[0][2])      # baseline model
        y = np.add(y, gaussian_peak(x, para[1]))
        y = np.add(y, gaussian_peak(x, para[2]))
        return y

    def gaussian_2peak(x, para): # para: a, b, c, d; A1, A2, B1, C1, B2, C2
        y = np.add(np.add(np.multiply(para[0], np.square(np.subtract(x, para[3]))), np.multiply(para[1], np.subtract(x, para[3]))), para[2])       # baseline model
        p1 = [para[4]]+para[6:8]
        p2 = [para[5]]+para[8:10]
        y = np.add(y, gaussian_peak(x, p1))
        y = np.add(y, gaussian_peak(x, p2))
        return y

    def objective_globle(x, para):
        x0 = x[0:200] # pull out each x if they were stacked
        x1 = x[200:400]
        x2 = x[400:600]
        x3 = x[600:800]
        x4 = x[800:1000]
        p0 = para[0:10]
        p1 = para[10:16] + para[6:10]
        p2 = para[16:22] + para[6:10]
        p3 = para[22:28] + para[6:10]
        p4 = para[28:34] + para[6:10]
        y0 = gaussian_2peak(x0, p0)
        y1 = gaussian_2peak(x1, p1)
        y2 = gaussian_2peak(x2, p2)
        y3 = gaussian_2peak(x3, p3)
        y4 = gaussian_2peak(x4, p4)
        y = list(y0) + list(y1) + list(y2) + list(y3) + list(y4)
        return y


    def fit_global(x0, x1, x2, x3, x4, y0, y1, y2, y3, y4, para_guess, bounds, option):
        x_value = list(x0) + list(x1) + list(x2) + list(x3) + list(x4)
        y_value = list(y0) + list(y1) + list(y2) + list(y3) + list(y4)
        para_g = para_guess + para_guess[0:6] + para_guess[0:6] + para_guess[0:6] + para_guess[0:6]
        bounds_g = {'ub':[], 'lb':[]}
        bounds_g['ub'] = bounds['ub'] + bounds['ub'][0:6] + bounds['ub'][0:6] + bounds['ub'][0:6] + bounds['ub'][0:6]
        bounds_g['lb'] = bounds['lb'] + bounds['lb'][0:6] + bounds['lb'][0:6] + bounds['lb'][0:6] + bounds['lb'][0:6]

        y = objective_globle(x_value, para_g)
        pyplot.figure(500)
        pyplot.plot(y)
        # pyplot.show()

        results_global = pyjcfit(objective_globle, x_value, y_value, para_g, bounds_g, option)
        results = results_global
        yfit = objective_globle(x_value, results['para'])
        residual = np.subtract(y_value, yfit)
        results['yfit'] = yfit
        results['residual'] = residual

        pyplot.figure()
        pyplot.subplot(2, 1, 1)
        x_new = range(len(y_value))
        pyplot.scatter(x_new,y_value)
        pyplot.plot(x_new, yfit, '-', color='red')
        pyplot.plot(x_new, residual, '-', color='orange')
        pyplot.title('pyjcfit: ')
        # pyplot.show()

        pyplot.subplot(2, 1, 2)
        pyplot.plot(results['error_hist'], '-', color='red')
        pyplot.title('error history')

        return results


    # or simulate a set of data
    n = 200
    para_true = ((-0.00002, 0.001, 0.5, 0), (1, 98, 15), (1, 103, 8)) # (a, b, c, d),  (a1, b1, c1), (a2, b2, c2)
    x_value = np.arange(0, n, 1)
    # y_value = list(objective(x_value, para_true))
    x_valuen = x_value + np.random.normal(0, 0.1, len(x_value))
    y_value = list(np.add(data_simulation(x_valuen, para_true), np.random.normal(0, 0.01, len(x_value)))) #SNR about sum of paraTrue(A)/noise

    #  simulate t = 0-4 with [C] = exp(-kt) and k = 0.30
    x0 = x_value + np.random.normal(0, 0.1, len(x_value))
    x1 = x_value + np.random.normal(0, 0.1, len(x_value))
    x2 = x_value + np.random.normal(0, 0.1, len(x_value))
    x3 = x_value + np.random.normal(0, 0.1, len(x_value))
    x4 = x_value + np.random.normal(0, 0.1, len(x_value))
    para0 = ((-0.000021, 0.0011, 0.51, 0), (1, 95, 15), (0, 105, 8)) # changing A1 and A2 and slightly
    para1 = ((-0.000022, 0.0012, 0.52, 10), (0.74, 95, 15), (0.26, 105, 8))  # changing A1 and A2 and slightly
    para2 = ((-0.000023, 0.0013, 0.53, 20), (0.55, 95, 15), (0.45, 105, 8))  # changing A1 and A2 and slightly
    para3 = ((-0.000024, 0.0014, 0.54, 30), (0.41, 95, 15), (0.59, 105, 8))  # changing A1 and A2 and slightly
    para4 = ((-0.000025, 0.0015, 0.55, 40), (0.30, 95, 15), (0.70, 105, 8))  # changing A1 and A2 and slightly
    y0 = list(np.add(data_simulation(x0, para0), np.random.normal(0, 0.01, len(x_value))))
    y1 = list(np.add(data_simulation(x1, para1), np.random.normal(0, 0.01, len(x_value))))
    y2 = list(np.add(data_simulation(x2, para2), np.random.normal(0, 0.01, len(x_value))))
    y3 = list(np.add(data_simulation(x3, para3), np.random.normal(0, 0.01, len(x_value))))
    y4 = list(np.add(data_simulation(x4, para4), np.random.normal(0, 0.01, len(x_value))))

    pyplot.figure()
    pyplot.plot(x0,y0)
    pyplot.plot(x1, y1)
    pyplot.plot(x2, y2)
    pyplot.plot(x3, y3)
    pyplot.plot(x4, y4)
    # pyplot.show()



# ------JCFIT starts here initial guess of parameters and bounds
    print('\n---- JCFit starts ---- ')
    para_guess = [-0.00002, 0.001, 0.5, 0,1, 1, 100, 10, 100, 10] # initial guessed parameters for the objective function
    bounds = {'ub': (1, 1, 100, 100, 2, 2, 200, 30, 200, 30), 'lb': (-1, -1, -100, -100, 0, 0, 1, 5,  1, 5)}  # upper bounds (ub) and lower bounds (lb) ordered the same as parameters
    # check if bounds and para_guess are well-ordered.
    d1 = np.subtract(bounds['ub'], para_guess)
    d2 = np.subtract(para_guess, bounds['lb'])
    d3 = np.multiply(d1, d2)
    for i in range(len(d3)):
        if d3[i] < 0:
            print(" para_guess[%s] and its bounds out of order." % i)

    # set searching options
    option = {'maxiteration': 500, 'precision': 1E-6, 'exp_step': 0.5, 'convgtest': 0, 'score_func': 'L2', 'show_error_surface': False}
    # maxiteration is the maximum searching iteration.
    # precision defines the significant figures. It is the smallest numerical search step of each paramter. e.g. paraguess of previous iteration = 10 and precision = 0.01, then searching step is 0.1 for this iteration and this parameter, i.e. precision = 0.01 is 2 sig fig.
    # exp_step, searching step size +-para*precision*(2^exp_step)^n where n is 1, 2, 3,...
    # convgtest is the minimum difference of sum(residual**2) between searching steps to judge convergence.

    #  --------- fitting starts
    start_time = time.time()
    fit_results = fit_global(x0, x1, x2, x3, x4, y0, y1, y2, y3, y4, para_guess, bounds, option)
    # fit_results = pyjcfit(objective, x_value, y_value, para_guess, bounds, option)
    finish_time = time.time()


    print('\n---- elapsed time for pyjcfit = %.9f seconds ----' % (finish_time - start_time))

    para = fit_results['para']
    para_hist = fit_results['para_hist']
    error_hist = fit_results['error_hist']
    yfit = fit_results['yfit']
    residual = fit_results['residual']
    # residual = np.subtract(y_value, objective(x_value, para))
    goodness_of_fit = fit_results['gof']

    # ------- fitting ends-----
    #
    # # ------ plot results
    print('\nfitted para = ', end = '')
    for i in para:
        print('%.5f' %i, end = ' ')
    print('\nTrue para   = ', end='')
    print('para_true:  ', para_true)
    print('\ngoodness of fitting: ', goodness_of_fit)

    A1 = [para[4], para[14], para[20], para[26], para[32]]
    A2 = [para[5], para[15], para[21], para[27], para[33]]


    def exp_decay(x, b, a, k):  # for scipy curve_fit
        return b + a * np.exp(x * (-k))
    para1, _ = curve_fit(exp_decay, range(5), A1, [1, 1, 1])
    para2, _ = curve_fit(exp_decay, range(5), A2, [1, 1, 1])
    x_new = np.arange(0, 5, 0.1)
    y_new1 = exp_decay(x_new, para1[0], para1[1], para1[2])
    y_new2 = exp_decay(x_new, para2[0], para2[1], para2[2])
    print('P1: ', para1, '  P2: ', para2)

    pyplot.figure()
    pyplot.scatter(range(5), A1)
    pyplot.scatter(range(5), A2)
    pyplot.plot(x_new, y_new1)
    pyplot.plot(x_new, y_new2)
    pyplot.title(['P1: ', para1, ' P2: ', para2])


    pyplot.show()

# end of example Jixin Chen @Ohio University, Athens, Ohio 11/2023

