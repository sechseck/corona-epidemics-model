#!/usr/bin/env python3

import sys
import argparse
import pprint

from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import pandas as pd
from numpy.linalg import norm
from numpy import asarray,hstack
import matplotlib.pyplot as plt


def loadData(region):
    """Load case data for region from CSV file

    :param region: Name of the region for which to load data

    :returns: Case data for region as pandas DataFrame
    """
    filename = r'./case_numbers.csv'
    data = pd.read_csv(filename, delimiter=';', comment="#")
    C = data[data['region']==region]['cases'].to_numpy()

    # If no cases are found, throw an error
    if len(C) == 0:
        print("Found 0 datapoints for region {}, terminating".format(region))
        raise SystemExit(3)
    return C


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

    :param t: TODO

    :param SIR: TODO

    :param beta: Model transition probability

    :param gamme: Model transisiton probability

    :returns: TODO
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
    """Solve differnetial equation TODO

    :param bgS: TODO

    :param IRC: TODO

    :returns: TODO
    """
    C = IRC[2:]
    tmax = len(C)
    (t, S, I, R) = solveode((bgS[2], IRC[0], IRC[1]), (bgS[0], bgS[1], tmax))
    deviation = norm(C - (I + R))
    return deviation


def solveode(SIR,bgt):
    """Solve the sir equations for the time-dependent S,I,R values

    :param C: Observed case numbers

    :param bgt: TODO

    :returns: TODO
    """
    S0 = SIR[0]
    I0 = SIR[1]
    R0 = SIR[2]
    beta = bgt[0]
    gamma = bgt[1]
    tmax = bgt[2]

    # Solve the SIR ode for S,I,R, given a beta and gamma
    sir_sol = solve_ivp(sir_ode,(0,tmax),(S0,I0,R0),args=(beta,gamma), \
                        t_eval=range(tmax))
    t =  sir_sol['t']
    S = sir_sol['y'][0]
    I = sir_sol['y'][1]
    R = sir_sol['y'][2]
    return (t, S, I, R)


def plot_results(results, C):
    """Create a plot visualizing curves for S, I and R

    :param results: Parameter estimation results

    :param C: Observed case numbers
    """
    fig, ax = plt.subplots(1,3)
    tmax = 2*len(C)
    (t, S, I, R) = solveode((results['S'],results['I0'],results['R0']), \
                            (results['beta'],results['gamma'],tmax))

    # Plot data and set axis subplot titles
    ax[0].plot(t,S)
    ax[0].set_title('S')
    ax[1].plot(t,I)
    ax[1].set_title('I')
    ax[2].plot(t,R)
    ax[2].set_title('R')
    plt.show()


if __name__ == "__main__":

    # Script only tested with Python > 3.7
    if(not sys.version_info >= (3, 7)):
        print("This script requires at least Python 3.7, terminating.")
        sys.exit(1)

    # Parse commmand line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", metavar="REGION", default="china", \
                        dest="region",
                        help="specify region for which to solve model")
    args = parser.parse_args()

    # Estimate parameters for SIR model and display results
    region = args.region
    C = loadData(region)
    results = parmest(C)
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(results)
    plot_results(results, C)
