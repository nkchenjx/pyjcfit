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

def pyjcfit(f, xdata, ydata, paraguess, bounds = {}, option = {'maxiteration': 50, 'precision': 0.00001, 'exp_step': 0.5, 'convgtest': 1E-100}):
    para = paraguess.copy()
    if len(bounds) == 0:
        ub =  [1E10] * len(para)
        lb = [-1E10] * len(para)
        bounds = {'ub': ub, 'lb': lb}
    para_hist = [para.copy()] * (option['maxiteration'] + 1)
    error_hist = [0] * option['maxiteration']
    errorlast = 0
    for iteration in range(option['maxiteration']):
        if (iteration + 1) % 100 == 0:
            print('.')
        else:
            print('.', end='')

        for i in range(len(para)):
            p = para[i]
            precision = option['precision'] * (abs(para[i]) + option['precision'])
            lb = bounds['lb'][i]
            ub = bounds['ub'][i]
            ll = p - lb
            nl = np.arange(int(np.floor(np.log2(ll / precision + 1))), 1, -option['exp_step'])
            ul = ub - p
            nu = np.arange(1, int(np.floor(np.log2(ul / precision + 1))), option['exp_step'])
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
        para_hist[iteration + 1] = para.copy()
        error_hist[iteration] = error[indmin].copy()
        # convergence test
        if abs(error[indmin] - errorlast) <= option['convgtest']:
            print('\n convergence reached')
            break
        errorlast = error[indmin]
    # print(para_hist)
    gof = pyjcfit_goodness(f, xdata, ydata, para, bounds, option)
    return {'para': para, 'para_hist': para_hist, 'error_hist': error_hist, 'gof': gof}


def pyjcfit_goodness(f, xdata, ydata, para, bounds, option):
    yfit = f(xdata, para)
    n = len(xdata)
    residual = np.subtract(ydata, yfit)
    residual_sq = np.square(residual, residual)
    a = np.divide(residual_sq, yfit)
    chisq = sum(filter(lambda i: not(np.isinf(i)), a))

    ymean =np.mean(yfit)
    sumy = sum(np.square(np.subtract(ydata, ymean)))
    sumr = sum(residual_sq)
    rsq = 1-sumr/sumy
    sigma = np.sqrt(sumr/n)

    # find the 95% confidence region of the parameters
    para95_ub = [0]*len(para)
    para95_lb = [0]*len(para)

    for i in range(len(para)):
        precision = option['precision'] * (abs(para[i]) + option['precision'])
        parau = para.copy()
        p = para[i]
        ub = bounds['ub'][i]
        ul = ub-p
        nu = np.arange(1, int(np.floor(np.log2(ul / precision + 1))), option['exp_step'])
        psu = np.append(p, np.add(p, np.multiply(np.exp2(nu), precision)))
        sigmau = [0]*len(psu)
        for j in range(len(psu)):
            parau[i] = psu[j]
            yfitu = f(xdata, parau)
            residual = np.subtract(ydata, yfitu)
            sigmau[j] = np.sqrt(np.dot(residual, residual)/n)
        a = abs(np.subtract(sigmau , sigma))
        b = abs(a-2*sigma/np.sqrt(n-len(para)))
        indminu = np.array(b).argmin()
        para95_ub[i] = psu[indminu]

        paral = para.copy()
        lb = bounds['lb'][i]
        ll = p - lb
        nl = np.arange(int(np.floor(np.log2(ll / precision + 1))), 1, -option['exp_step'])
        psl = np.append(np.subtract(p, np.multiply(np.exp2(nl), precision)), p)
        sigmal = [0] * len(psl)
        for j in range(len(psl)):
            paral[i] = psl[j]
            yfitl = f(xdata, paral)
            residual = np.subtract(ydata, yfitl)
            sigmal[j] = np.sqrt(np.dot(residual, residual)/n)
        a = abs(sigmal - sigma)
        b = abs(a-2*sigma/np.sqrt(n-len(para)))
        indminl = np.array(b).argmin()
        para95_lb[i] = psl[indminl]

        # # show the error vs parameter curve
        # pyplot.plot(psu, sigmau)
        # pyplot.plot(psl, sigmal)
        # pyplot.scatter(psu[indminu], sigmau[indminu])
        # pyplot.scatter(psl[indminl], sigmal[indminl])
        # pyplot.show()


    goodness_of_fit = {'rsq': rsq, 'chisq': chisq, 'sigma': sigma, 'para95_lb': para95_lb, 'para95_ub': para95_ub}
    return goodness_of_fit


# ----------- An example objective function and main function  ---------------:
def objective(x, para): #must have all parameters in the list para
    return para[0] * x + para[1]


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
        print("bounds unusual, double check bounds and paraguess...")
    # set searching options
    option = {'maxiteration': 50, 'precision': 0.0001, 'exp_step': 0.5, 'convgtest': 1E-100}
    # maxiteration is the maximum searching iteration.
    # precision defines the significant figures. It is the smallest numerical search step of each paramter. e.g. paraguess of previous iteration = 10 and precision = 0.01, then searching step is 0.1 for this iteration and this parameter, i.e. precision = 0.01 is 2 sig fig.
    # exp_step, searching step size +-para*precision*(2^exp_step)^n where n is 1, 2, 3,...
    # convgtest is the minimum difference of sum(residual**2) between searching steps to judge convergence.

    #  --------- fitting starts
    fitresults = pyjcfit(objective, x_value, y_value, paraguess, bounds, option)
    para = fitresults['para']
    para_hist = fitresults['para_hist']
    error_hist = fitresults['error_hist']
    residual = np.subtract(y_value, objective(x_value, para))
    goodness_of_fit = fitresults['gof']
    #  -------_ fitting ends

    #------ plot results
    print('\nfitted para = ', para)
    print('goodness of fitting: ', goodness_of_fit)

    x_new = np.arange(0, 20, 0.1)
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
