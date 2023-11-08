# Coded by Jixin Chen @ Ohio University, Department of Chemistry and Biochemistry
# First coded in MATLAB on 2016/04/05
# Converted to python 3.12 2023/11/07

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


# from scipy.optimize import curve_fit

def pyjcfit(f, xdata, ydata, paraguess, bounds, option):
    para = paraguess
    para_hist = [paraguess] * (option['maxiteration'] + 1)
    error_hist = [0] * option['maxiteration']
    errorlast = 0
    for iteration in range(option['maxiteration']):
        if (iteration + 1) % 100 == 0:
            print('\n')
        else:
            print('.', end='')

        for i in range(len(para)):
            p = para[i]
            precision = option['precision']
            lb = bounds['lb'][i]
            ub = bounds['ub'][i]
            ll = p - lb
            nl = range(int(np.floor(np.log2(ll / precision + 1) * 2)), 2, -1)
            nl = np.divide(list(nl), 2)
            ul = ub - p
            nu = range(2, int(np.floor(np.log2(ul / precision + 1) * 2)), 1)
            nu = np.divide(list(nu), 2)
            ps = np.append(lb, np.subtract(p, np.multiply(np.exp2(nl), precision)))
            ps = np.append(ps, np.add(p, np.multiply(np.exp2(nu), precision)))
            ps = np.append(ps, ub)
            error = [0] * len(ps)
            for j in range(len(ps)):
                para[i] = ps[j]
                residual = ydata - f(xdata, para)
                error[j] = np.dot(residual, residual)
            indmin = np.array(error).argmin()
            para[i] = ps[indmin]
        para_hist[iteration + 1] = para
        error_hist[iteration] = error[indmin]
        # convergence test
        if abs(error[indmin] - errorlast) <= option['convgtest']:
            print('\n convergence reached')
            break
        errorlast = error[indmin]
    # print(para_hist)
    fe = pyjcfit_finderror(f, xdata, ydata, para, bounds, option)
    return {'para': para, 'para_hist': para_hist, 'error_hist': error_hist}


def pyjcfit_finderror(f, xdata, ydata, para, bounds, option):
    yfit = f(xdata, para)
    residual = np.divide(ydata, yfit)

    return


# ----------- An example objective function and main function  ---------------
def objective(x, para):
    return para[0] * x + para[1]


# def objectivec(x, a, b):  # for scipy curve_fit
#     return a * x + b

if __name__ == '__main__':
    # load input variables

    data = pd.read_csv("data.csv")
    # print(data)
    x_value = data.x
    y_value = data.y
    #   print(y_value)

    # initial guess of parameters and bounds
    paraguess = [1, 1]  # parameters for the objective function
    bounds = {'ub': [1000, 1000], 'lb': [-1000, -1000]}  # bounds ordered the same as paramters
    # check if bounds and paraguess are well ordered.
    d1 = np.subtract(bounds['ub'], paraguess)
    d2 = np.subtract(paraguess, bounds['lb'])
    d3 = np.multiply(d1, d2)
    d4 = [i for i in d3 if i < 0]
    if d4:
        print("bounds unusual")
    # set searching options
    option = {'maxiteration': 50, 'precision': 0.0001, 'convgtest': 1E-100}
    # maxiteration is the maximum searching iteration.
    # precision is the smallest numerical search step.
    # convgtest is the minimum difference of sum(residual**2) between searching steps to judge convergence.

    #  --------- fitting starts
    # para, _ = curve_fit(objectivec, x_value, y_value)
    fitresults = pyjcfit(objective, x_value, y_value, paraguess, bounds, option)
    para = fitresults['para']
    para_hist = fitresults['para_hist']
    error_hist = fitresults['error_hist']
    residual = np.divide(y_value, objective(x_value, para))
    #  -------_ fitting ends

    # plot results
    x_new = np.arange(0, 20, 0.1)
    print('\nfitted para = ', para)
    y_new = objective(x_new, para)

    pyplot.subplot(1, 2, 1)
    pyplot.scatter(x_value, y_value)
    pyplot.plot(x_new, y_new, '-', color='red')
    pyplot.plot(x_value, residual, '--', color='orange')
    pyplot.title(str('y =%.5f * x + (%.5f)' % (para[0], para[1])))
    # pyplot.show()

    pyplot.subplot(1, 2, 2)
    pyplot.plot(error_hist, '-', color='red')
    pyplot.title('error history')
    pyplot.show()
# end of example
