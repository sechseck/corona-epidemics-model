#!/usr/bin/env python3

import sys
import argparse
import pprint

from scipy.integrate import solve_ivp
from scipy.optimize import minimize,curve_fit
import pandas as pd
from numpy.linalg import norm
from numpy import asarray,hstack,exp,log,arange,concatenate
import matplotlib.pyplot as plt
from matplotlib import rcParams


def loadData(country):
    """Load case data for country from CSV file

    :param country: Name of the region for which to load data

    :returns: Case data for region as pandas DataFrame
    """
    filename = r'./case_numbers.csv'
    data = pd.read_csv(filename, delimiter=';', comment="#")
    C = data[data['country'] == country]['cases'].to_numpy()

    # If no cases are found, throw an error
    if len(C) == 0:
        print("Found 0 datapoints for region {}, terminating".format(country))
        raise SystemExit(3)
    return C

def stubData(C):
    """create simulated data based on exponential growth for the first days if first case count is > 1, i.e. first cases have been unobserved

    :param C: the observed case numbers
    :return: stubC: the simulated case numbers, (a,b,c): the parameters of the exponential function
    """

    trainingPeriodIdxMax = 5 #the first trainingPeriodIdxMax case numbers will be used to estimate the function (or less if not available)
    first_sim_case_count = 80
    trainingPeriod = range(trainingPeriodIdxMax if len(C) > trainingPeriodIdxMax else len(C))

    #fit the function to the data
    popt, pcov = curve_fit(funcExponential, trainingPeriod, C[trainingPeriod])
    (a,b,c) = popt
    #calculate the (negative) simulated date number where case count = first_sim_case_count
    firstSimDateNumber = int(invFuncExponential(first_sim_case_count,a,b,c))
    simDateNumbers = arange(firstSimDateNumber,0)

    result = funcExponential(simDateNumbers, *popt).astype(int)
    parms = {'a':a,'b':b,'c':c}


    return {'result': result, 'parms': parms}

def funcExponential(x, a, b, c):
    """ the generic exponential function used in the curve fitting routine
    :param x: list of x (datenumber) values for function evaluation
    :param a: value of function at x == offset
    :param b: growth factor
    :param c: x offset
    :return:
    """
    return a * exp(b *(x-c))

def invFuncExponential(y, a, b, c):
    """the inverse of the exponential function
    :param y: list of y (case count) values for function evaluation
    :param a:
    :param b:
    :param c:
    :return:
    """
    return 1 / b * log(y / a) + c


def parmest(C):
    """Estimate model parameters based on data

    This function fits the observation data to the model.

    :param C: The time series of observed cases

    :returns: Parameters estimates for 'beta', 'gamma', 'S',
        'I0', 'R0' in a dict
    """
    firstCaseCount = C[0]
    [betaGuess,gammaGuess,SGuess,I0,R0] = iniguess(firstCaseCount)

    # Mangle I0, R0 and array C into 1-d ndarray
    parmArray = hstack((asarray([I0,R0]),C))

    optOptions = {'maxiter': 20000, 'disp': False}

    #Find beta, gamma and S that fit the data C best given I0,R0
    optRes = minimize(optSolveOde,(betaGuess,gammaGuess,SGuess), \
                      method = 'Nelder-Mead', args=parmArray, \
                      options=optOptions)
    if optRes.success:
        [beta, gamma, S] = optRes.x
    else:
        print('Optimzation did not converge, terminating')
        raise SystemExit(2)

    return {'beta':beta, 'gamma':gamma, 'S':S, 'I0': I0, 'R0': R0}


def iniguess(firstCaseCount):
    """Guess initial parameters for parameter estimation
    :param firstCaseCount: Number of infected at time t0
    :returns: Initial parameters guess for 'beta', 'gamma', 'S',
        'I0', 'R0' in a dict
    """
    betaGuess = 10
    gammaGuess = 10

    # Number of susceptible people at t0. normalized to one,
    # but using a big number for numerical stability
    SGuess = 1e9-firstCaseCount

    # Number of infected at t0
    I0 = firstCaseCount

    # Number of recovered at t0
    R0 = 0

    return (betaGuess,gammaGuess,SGuess,I0,R0)


def sir_ode(t,SIR,beta,gamma):
    """The SIR  differential equation
    :param t: the time variable, not used
    :param SIR: current values of S,I,R
    :param beta: Model transition probability
    :param gamme: Model transisiton probability
    :returns: solutions for dS,dI,dR as functions of time
    """
    S = SIR[0]
    I = SIR[1]
    R = SIR[2]
    N = S + I + R

    # Define ODEs
    dS = -beta*I*S/N
    dI = beta*S*I/N - gamma*I
    dR = gamma*I
    return [dS,dI,dR]


def optSolveOde(bgS,IRC):
    """Solve differential equation and compute deviation from observed values

    :param bgS: beta,gamma,S are variables for the differential equation

    :param IRC: I,R and measurements C are treated as constant

    :returns: the euclidean norm of the differences between model estimates and observed values
    """
    C = IRC[2:]
    tmax = len(C)
    (t, S, I, R) = solveode((bgS[2], IRC[0], IRC[1]), (bgS[0], bgS[1], tmax))
    deviation = norm(C - (I + R))
    return deviation


def solveode(SIR,bgt):
    """Solve the differential equations

    :param SIR: time-dependent S,I,R values

    :param bgt: beta, gamma and maximum date tmax

    :returns: time index t, S,I,R as time-dependent function values
    """
    S0 = SIR[0]
    I0 = SIR[1]
    R0 = SIR[2]
    beta = bgt[0]
    gamma = bgt[1]
    tmax = bgt[2]

    # Solve the SIR ode for S,I,R, given a beta and gamma
    sir_sol = solve_ivp(sir_ode,(0,tmax),(S0,I0,R0),args=(beta,gamma),
                        t_eval=range(tmax))
    t =  sir_sol['t']
    S = sir_sol['y'][0]
    I = sir_sol['y'][1]
    R = sir_sol['y'][2]
    return (t, S, I, R)


def plot_SIR(t, S, I, R):
    """Create a plot visualizing curves for S, I and R
    :param t: date numbers for S,I,R
    :param S: number of susceptible individuals as a function of t
    :param I: number of infected individuals as a function of t
    :param R: number of recovered (healthy or dead) as a function of t
    :return: nothing, creates a plot
    """

    fig, ax = plt.subplots(1,3)

    # Plot data and set axis subplot titles
    ax[0].plot(t,S)
    ax[0].set_title('S(usceptible)')
    ax[0].set_xlabel('day count')
    ax[0].set_ylabel('number of individuals')
    ax[1].plot(t,I)
    ax[1].set_title('I(nfected)')
    ax[1].set_xlabel('day count')
    ax[1].set_ylabel('number of individuals')
    ax[2].plot(t,R)
    ax[2].set_title('R(ecovered or Dead)')
    ax[2].set_xlabel('day count')
    ax[2].set_ylabel('number of individuals')

def plot_goodness_of_fit(t, S, I, R, tC, C, region):
    """Create a plot visualizing curves for S, I and R
    :param t: date numbers for S,I,R
    :param S: number of susceptible individuals as a function of t, not needed here
    :param I: number of infected individuals as a function of t
    :param R: number of recovered (healthy or dead) as a function of t
    :param tC: date numbers for C
    :param C: observed case numbers
    :return: nothing, creates a plot
    """
    fig, ax = plt.subplots()

    markerSize = 2.6*rcParams['lines.markersize'] ** 2
    ax.scatter(tC,C,s=markerSize,c='red')
    ax.plot(t,I+R)

    ax.set_title('predicted and observed case numbers for\n' + region)
    ax.set_xlabel('day count')
    ax.set_ylabel('case number')


if __name__ == "__main__":

    # Script only tested with Python > 3.7
    if(not sys.version_info >= (3, 7)):
        print("This script requires at least Python 3.7, terminating.")
        sys.exit(1)

    # Parse commmand line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", metavar="COUNTRY", default="CHINA",
                        dest="country",
                        help="specify country for which to solve model")
    args = parser.parse_args()


    #configuration section
    minFirstCases = 10 # we will create simulated stub data if the first case count is below this figure #TODO: set this in relation to the current figure, i.e. a ratio


    # load data
    country = args.country
    C = loadData(country)

    #if the first value of C is too high, create simulated previous values
    if C[0] > minFirstCases:
        need2Stub = True
        stubResults = stubData(C)
        stubC = stubResults['result']
        finalC = concatenate((stubC, C))
    else:
        need2Stub = False
        finalC = C

    # calcuate results
    results = parmest(finalC) #estimate the model
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(results)

    # analyze and graph results
    #solve for S,I,R given the estimated parameters Expectation(S0),beta,gamma
    tmax = 2*len(C)
    (t, S, I, R) = solveode((results['S'],results['I0'],results['R0']),
                            (results['beta'],results['gamma'],tmax))
    plot_SIR(t, S, I, R)
    plot_goodness_of_fit(t, S, I, R, range(len(C)), C,country)
    plt.show()
