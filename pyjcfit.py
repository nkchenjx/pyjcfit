# Use a non-linear least square random searching algorithm to fit data.
# Coded by Jixin Chen @ Ohio University, Department of Chemistry and Biochemistry
# First coded in MATLAB on 2016/04/05
# Converted to python 3.12 2023/11/07. Developed and tested on software versions Python 3.12, Windows 11, Anaconda 3 11/2023, in PyCharm 2023.2.
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
#
#
# def pyjcfit(f, xdata, ydata, para_guess, bounds = {}, option = {'maxiteration': 50, 'precision': 0.00001, 'exp_step': 0.5, 'convgtest': 1E-100}):
# Use non-linear least square random searching algorithm to fit a function f to data.
# assume ydata = f(xdata, parameters)
# :param f: objective function input xdata and para, give yfit. can be complicated combination of functions of
#        multi-dimension or global data structure then vectorize x and y.
# :param xdata: independent variable, array-like. For high dimension or multiple curves, stack and vectorize
# :param ydata: dependent variable, array-like. If results are nD or a series of data, stack and vectorize
# :param para_guess: initial guess of the parameters of the model f. All parameters in a vector array, e.g. list.
# :param bounds: bounds of the parameters. can give a very large bound but the narrower the range the faster the fitting
# :param option: maxiteration, maximum number of iteration, necessary for most fitting projects to control time.
#        precision, significant figures of parameters. ex_step, searching spacing exponentially distributed,
#        e.g. 0.5 for a parameter 1.3 and precision 0.1 meaning the searching steps will be 1.3+2^0.5*0.1, 1.3+2^1*0.1, ...
# :return: {'para': para, 'para_hist': para_hist, 'error_hist': error_hist, 'gof': gof}
#          para is the fitted parameter
#          para_hist is the history of the parameter of the searching iterations.
#          error_hist is the square of residual over the iterations.
#          gof, goodness of fitting contains R-square, chi-square, sigma of the residual,
#               and bounds of the parameters with 95% confidence
#
# References
# ----------
# [1] https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.6b05697
#     Jixin Chen, Joseph R Pyle, Kurt Waldo Sy Piecco, Anatoly B Kolomeisky, Christy F Landes,
#     A Two-Step Method for smFRET Data Analysis, J. Phys. Chem. B, 2016, 120 (29), pp 7128–7132
# [2] Juvinch R. Vicente, Ali Rafiei Miandashti, Kurt Sy Piecco, Joseph R. Pyle, Martin E. Kordesch, Jixin Chen*,
#     Single-Particle Organolead Halide Perovskite Photoluminescence as a Probe for Surface Reaction Kinetics.
#     ACS Applied Matierals & Interfaces, 2019, 11(19), 18034-18043.
#
# # ----------- An example objective function and main function  ---------------:
# if __name__ == '__main__':
#     import numpy as np
#     import pandas as pd
#     from matplotlib import pyplot
#     from pyjcfit import pyjcfit
#     import time
#
#     def objective(x, para):
#         return para[0] * x + para[1]
#
#
#     # load input variables
#     data = pd.read_csv("data.csv")
#     # print(data)
#     x_value = data.x
#     y_value = data.y
#     #   print(y_value)
#
#     # initial guess of parameters and bounds
#     para_guess = [1, 1]  # initial guessed parameters for the objective function
#     bounds = {'ub': [1000, 1000],
#               'lb': [-1000, -1000]}  # upper bounds (ub) and lower bounds (lb) ordered the same as parameters
#     # check if bounds and para_guess are well-ordered.
#     d1 = np.subtract(bounds['ub'], para_guess)
#     d2 = np.subtract(para_guess, bounds['lb'])
#     d3 = np.multiply(d1, d2)
#     for i in range(len(d3)):
#         if d3[i] < 0:
#             print(" para_guess[%s] and its bounds out of order." % i)
#
#     # set searching options
#     option = {'maxiteration': 50, 'precision': 0.001, 'exp_step': 0.5, 'convgtest': 1E-100}
#     # maxiteration is the maximum searching iteration.
#     # precision defines the significant figures. It is the smallest numerical search step of each paramter. e.g. paraguess of previous iteration = 10 and precision = 0.01, then searching step is 0.1 for this iteration and this parameter, i.e. precision = 0.01 is 2 sig fig.
#     # exp_step, searching step size +-para*precision*(2^exp_step)^n where n is 1, 2, 3,...
#     # convgtest is the minimum difference of sum(residual**2) between searching steps to judge convergence.
#
#     #  --------- fitting starts
#     start_time = time.time()
#     fit_results = pyjcfit(objective, x_value, y_value, para_guess, bounds, option)
#     print('\n---- elapsed time for pyjcfit = %f seconds ----' % (time.time() - start_time))
#     para = fit_results['para']
#     para_hist = fit_results['para_hist']
#     error_hist = fit_results['error_hist']
#     residual = np.subtract(y_value, objective(x_value, para))
#     goodness_of_fit = fit_results['gof']
#
#     # ------- fitting ends-----
#
#     # ------ plot results
#     print('\nfitted para = ', para)
#     print('goodness of fitting: ', goodness_of_fit)
#     print(fit_results['para_hist'])
#
#     x_new = np.arange(0, 20, 0.1)
#     y_new = objective(x_new, para)
#
#     pyplot.subplot(1, 2, 1)
#     pyplot.scatter(x_value, y_value)
#     pyplot.plot(x_new, y_new, '-', color='red')
#     pyplot.plot(x_value, residual, '--', color='orange')
#     pyplot.title(str('y =%.5f * x + (%.5f)' % (para[0], para[1])))
#     # pyplot.show()
#
#     pyplot.subplot(1, 2, 2)
#     pyplot.plot(error_hist, '-', color='red')
#     pyplot.title('error history')
#     pyplot.show()
#
# # end of example

import numpy as np
import pandas as pd
from matplotlib import pyplot
import random


# ---------The ENTRANCE-------------------:
def pyjcfit(f, xdata, ydata, para_guess, bounds={}, option={}):  # the main function to call
    """
    use non-linear least square random searching algorithm to fit a function f to data.
    assume ydata = f(xdata, parameters)
    :param f: objective function input xdata and para, give yfit. can be complicated combination of functions of
           multi-dimension or global data structure then vectorize x and y.
    :param xdata: independent variable, array-like. For high dimension or multiple curves, stack and vectorize
    :param ydata: dependent variable, array-like. If results are nD or a series of data, stack and vectorize
    :param para_guess: initial guess of the parameters of the model f. All parameters in a vector array, e.g. list.
    :param bounds: bounds of the parameters. can give a very large bound but the narrower the range the faster the fitting
    :param option: maxiteration, maximum number of iteration, necessary for most fitting projects to control time.
           precision, significant figures of parameters, e.g. 0.1 two sig. fig., 0.01, 3 sig. fig. et.. ex_step, searching spacing exponentially distributed,
           e.g. 0.5 for a parameter 1.3 and precision 0.1 meaning the searching steps will be 1.3+2^0.5*0.1, 1.3+2^1*0.1, ...
    :return: {'para': para, 'para_hist': para_hist, 'error_hist': error_hist, 'gof': gof}
             para is the fitted parameter
             para_hist is the history of the parameter of the searching iterations.
             error_hist is the square of residual over the iterations.
             gof, goodness of fitting contains R-square, chi-square, sigma of the residual,
                  and bounds of the parameters with 95% confidence

    An example is given in the _main_ function.

    References
    ----------
    [1] https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.6b05697
        Jixin Chen, Joseph R Pyle, Kurt Waldo Sy Piecco, Anatoly B Kolomeisky, Christy F Landes,
        A Two-Step Method for smFRET Data Analysis, J. Phys. Chem. B, 2016, 120 (29), pp 7128–7132
    [2] Juvinch R. Vicente, Ali Rafiei Miandashti, Kurt Sy Piecco, Joseph R. Pyle, Martin E. Kordesch, Jixin Chen*,
        Single-Particle Organolead Halide Perovskite Photoluminescence as a Probe for Surface Reaction Kinetics.
        ACS Applied Matierals & Interfaces, 2019, 11(19), 18034-18043.
    """
    if not bounds:
        ub = [1E20] * len(para_guess)
        lb = [-1E20] * len(para_guess)
        bounds = {'ub': ub, 'lb': lb}
    if not option:
        option = {'maxiteration': 100, 'precision': 1E-6, 'exp_step': 0.5, 'convgtest': 1E-100, 'score_func': 'L2',
                  'show_error_surface': False}
    else:
        if not 'score_func' in option:
            option.update({'score_func': 'L2'})
        if not 'maxiteration' in option:
            option.update({'maxiteration': 100})
        if not 'precision' in option:
            option.update({'precision': 0.00001})  # 5 sig. fig.
        if not 'exp_step' in option:
            option.update({'exp_step': 0.5})
        if not 'convgtest' in option:
            option.update({'convgtest': 1E-100})
        if not 'show_error_surface' in option:
            option.update({'show_error_surface': False})
    if option['score_func'] == 'L2':
        option['score_func'] = score_func_L2
    elif option['score_func'] == 'L1':
        option['score_func'] = score_func_L2
    else:
        print(
            'unknown scoring methods, please code it and change the option in pyjcfit. Right now will try least square.')
        option['score_func'] = score_func_L2
    # ----------done initialization of options----------

    # --------- call the kernel to fit the data---------
    results = pyjcfit_kernel(f, xdata, ydata, para_guess, bounds, option)
    # --------- calcualte the goodness of the fitting results---------

    results['gof'] = pyjcfit_goodness(f, xdata, ydata, results['para'], bounds, option)
    return results
    # 95% bounds in gof are 2 sigma uncertainty at local minima and do not necessarily cover the global minima


# --------- End of entrance function and return the project results --------------


# ------------------Kernel functions starts below------------------------------:
def score_func_L2(residual):
    return np.dot(residual, residual)  # sum the square of the residuals


def score_func_L1(residual):
    return np.abs(residual).sum()  # sum the absolute value of the residuals


def pyjcfit_kernel(f, xdata, ydata, para_guess, bounds, option):
    score_func = option['score_func']  # load the method name to calculate the scores
    para = para_guess.copy()
    para_hist = [None] * (option['maxiteration'] + 1)
    para_hist[0] = para.copy()
    error_hist = [None] * option['maxiteration']
    errorlast = 0
    error_counts = 0
    para_order = list(range(len(para)))
    
    for iteration in range(option['maxiteration']):
        if (iteration + 1) % 100 == 0:
            print('.')
        else:
            print('.', end='')
        random.shuffle(para_order)
        for i in para_order:
            p = para[i]
            precision = option['precision'] * (abs(para[i]) + option[
                'precision'])  # calculate the smallest step for the significant figure of the parameter
            lb = bounds['lb'][i]
            ub = bounds['ub'][i]
            ll = p - lb
            nl = np.arange(int(np.floor(np.log2(ll / precision + 1))), 1,
                           -option['exp_step'])  # calculate the exponential spacing sequence to the lower bound
            ul = ub - p
            nu = np.arange(1, int(np.floor(np.log2(ul / precision + 1))),
                           option['exp_step'])  # calculate the exponential spacing sequence to the upper bound
            ps = np.append(lb, np.subtract(p, np.multiply(np.exp2(nl), precision)))
            ps = np.append(ps, np.add(p, np.multiply(np.exp2(nu), precision)))
            ps = np.append(ps,
                           ub)  # done generate the searching series with exponentially distributed spaces for this parameter
            error = [0] * len(ps)
            for j in range(len(ps)):
                para[i] = ps[j]
                try:
                    residual = np.subtract(ydata, f(xdata, para))
                    error[j] = score_func(residual)
                except:
                    error[j] = np.inf  # in case of invalid calculations, set error inf.
                    error_counts += 1
            indmin = np.array(error).argmin()  # find where is the minimum error on this searching series then
            para[i] = ps[indmin]  # update the optimum
        para_hist[iteration + 1] = para.copy()
        error_hist[iteration] = error[indmin] / len(ydata)
        # convergence test
        if abs(error[indmin] - errorlast) <= option['convgtest']:
            print('\n convergence reached at # %i iteration' % iteration)
            break
        errorlast = error[indmin]
    # print(para_hist)
    if iteration == option['maxiteration'] - 1:
        print(
            'maximum iteration number reached, change it if needed and/or use the result para as initial guess to continue.')
    if error_counts:
        print('encountered ', error_counts, ' function calculation fails probably due to large bounds.')
    return {'para': para, 'para_hist': para_hist, 'error_hist': error_hist}


# ------- End of kernel functions----------------------------------------------------------------------

# ----------Check the goodness of the fitting results --------------------------
def pyjcfit_goodness(f, xdata, ydata, para, bounds, option):
    # score_func = option['score_func'] # will always use least-square to calculate the goodness of fitting
    yfit = f(xdata, para)
    n = len(xdata)
    residual = np.subtract(ydata, yfit)
    residual_sq = np.square(residual, residual)
    a = np.divide(residual_sq, yfit)
    chisq = sum(filter(lambda i: not (np.isinf(i)), a))

    ymean = np.mean(yfit)
    sumy = sum(np.square(np.subtract(ydata, ymean)))
    sumr = sum(residual_sq)
    rsq = 1 - sumr / sumy
    sigma = np.sqrt(sumr / n)

    # find the 95% confidence region of the parameters
    para95_ub = [0] * len(para)
    para95_lb = [0] * len(para)

    for i in range(len(para)):
        precision = option['precision'] * (abs(para[i]) + option['precision'])
        parau = para.copy()
        p = para[i]
        ub = bounds['ub'][i]
        ul = ub - p
        nu = np.arange(1, int(np.floor(np.log2(ul / precision + 1))), option['exp_step'])
        psu = np.append(p, np.add(p, np.multiply(np.exp2(nu), precision)))
        psu = np.append(psu, ub)
        sigmau = [0] * len(psu)
        for j in range(len(psu)):
            parau[i] = psu[j]
            try:
                yfitu = f(xdata, parau)
                residual = np.subtract(ydata, yfitu)
                sigmau[j] = np.sqrt(np.dot(residual, residual) / n)
            except:
                sigmau[j] = np.inf
        a = abs(np.subtract(sigmau, sigma))
        b = abs(a - 2 * sigma / np.sqrt(n - len(para)))
        indminu = np.array(b).argmin()
        para95_ub[i] = psu[indminu]

        paral = para.copy()
        lb = bounds['lb'][i]
        ll = p - lb
        nl = np.arange(int(np.floor(np.log2(ll / precision + 1))), 1, -option['exp_step'])
        psl = np.append(np.subtract(p, np.multiply(np.exp2(nl), precision)), p)
        psl = np.append(lb, psl)
        sigmal = [0] * len(psl)
        for j in range(len(psl)):
            paral[i] = psl[j]
            try:
                yfitl = f(xdata, paral)
                residual = np.subtract(ydata, yfitl)
                sigmal[j] = np.sqrt(np.dot(residual, residual) / n)
            except:
                sigma[j] = np.inf
        a = abs(sigmal - sigma)
        b = abs(a - 2 * sigma / np.sqrt(n - len(para)))
        indminl = np.array(b).argmin()
        para95_lb[i] = psl[indminl]

        # # show the error vs parameter curve
        try:
            if option['show_error_surface']:
                pyplot.figure()
                pyplot.plot(psu, sigmau)
                pyplot.plot(psl, sigmal)
                pyplot.scatter(psu[indminu], sigmau[indminu])
                pyplot.scatter(psl[indminl], sigmal[indminl])
                # pyplot.show()
        except Exception:
            print('show_error_surface is missing in option.')
    goodness_of_fit = {'rsq': rsq, 'chisq': chisq, 'sigma': sigma, 'para95_lb': para95_lb, 'para95_ub': para95_ub}
    # 95% bounds are 2 sigma uncertainty at local minima and do not necessarily cover the global minima
    return goodness_of_fit


# ----------- An example objective function and main function  ---------------:
if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot
    from pyjcfit import pyjcfit
    import time


    def objective(x, para):
        return np.add(np.multiply(para[0], x), para[1])


    # # load input variables
    # data = pd.read_csv("data.csv")
    # # print(data)
    # x_value = data.x
    # y_value = data.y
    # #   print(y_value)

    # or simulate a set of data
    n = 100
    para_true = (1.23456789, 9.87654321)
    x_value = np.arange(0, n, 10)
    x_valuen = x_value + np.random.normal(0, 0.00001, len(x_value))
    y_value = list(np.add(objective(x_valuen, para_true),
                          np.random.normal(0, 10, len(x_value))))  # SNR about sum of paraTrue(A)/noise

    # initial guess of parameters and bounds
    para_guess = [1, 1]  # initial guessed parameters for the objective function
    bounds = {'ub': (1000, 1000),
              'lb': (-1000, -1000)}  # upper bounds (ub) and lower bounds (lb) ordered the same as parameters
    # check if bounds and para_guess are well-ordered.
    d1 = np.subtract(bounds['ub'], para_guess)
    d2 = np.subtract(para_guess, bounds['lb'])
    d3 = np.multiply(d1, d2)
    for i in range(len(d3)):
        if d3[i] < 0:
            print(" para_guess[%s] and its bounds out of order." % i)

    # set searching options
    option = {'maxiteration': 100, 'precision': 1E-15, 'exp_step': 0.5, 'convgtest': 1E-100}
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

    pyplot.subplot(1, 2, 1)
    pyplot.scatter(x_value, y_value)
    pyplot.plot(x_new, y_new, '-', color='red')
    pyplot.plot(x_value, residual, '-', color='orange')
    pyplot.title(str('y =%.5f * x + (%.5f)' % (para[0], para[1])))
    # pyplot.show()

    pyplot.subplot(1, 2, 2)
    pyplot.plot(error_hist, '-', color='red')
    pyplot.title('error history')
    pyplot.show()

# end of example
