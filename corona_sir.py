from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import pandas as pd
from numpy.linalg import norm
from numpy import asarray,hstack
import matplotlib.pyplot as plt
import pprint

#load data from file
def loadData(country):
    filename = r'./fallzahlen.csv'
    data = pd.read_csv(filename, delimiter=';')
    C = data[data['country']==country]['cases'].to_numpy()
    return C

#parameter estimation from data, the routine that fits observation to model
#input params: C: the time series of observed cases
def parmest(C):
    firstCaseCount = C[0]
    [betaGuess,gammaGuess,SGuess,I0,R0] = iniguess(firstCaseCount)

    #mangle I0,R0 and array C into 1-d ndarray
    parmArray = hstack((asarray([I0,R0]),C))

    optOptions = {'maxiter': 20000, 'disp': False}
    optRes = minimize(optSolveOde,(betaGuess,gammaGuess,SGuess),method = 'Nelder-Mead',args=parmArray,options=optOptions) #find beta, gamma and S that fit the data C best given I0,R0
    if optRes.success:
        [beta, gamma, S] = optRes.x
    else:
        print ('Optimzation did not converge, terminating')
        raise SystemExit

    return {'beta':beta, 'gamma':gamma, 'S':S, 'I0': I0, 'R0': R0}

#guess initial parameters for parest parameter estimation
def iniguess(firstCaseCount):
    betaGuess = 1/0.00267103
    gammaGuess = 1/0.00267232
    SGuess = 1e9-firstCaseCount #number of susceptible people at t0. normalized to one, but using a big number for numerical stability
    I0 = firstCaseCount #number of infected at t0
    R0 = 0 #number of recovered at t0
    return (betaGuess,gammaGuess,SGuess,I0,R0)

#the sir differential equation
def sir_ode(t,SIR,beta,gamma):
    S = SIR[0]
    I = SIR[1]
    R = SIR[2]
    N = S + I + R
    # ODEs
    dS = -beta*I*S/N
    dI = beta*S*I/N-gamma*I
    dR = gamma*I
    return [dS,dI,dR]

def optSolveOde(bgS,IRC):
    C = IRC[2:]
    tmax = len(C)
    (t,S,I,R) = solveode((bgS[2],IRC[0],IRC[1]), (bgS[0],bgS[1],tmax))
    deviation = norm(C - (I + R))
    return deviation

#solve the sir equations for the time-dependent S,I,R values
# C: observed case numbers
def solveode(SIR,bgt):
    S0 = SIR[0]
    I0 = SIR[1]
    R0 = SIR[2]
    beta = bgt[0]
    gamma = bgt[1]
    tmax = bgt[2]
    sir_sol = solve_ivp(sir_ode,(0,tmax),(S0,I0,R0),args=(beta,gamma),t_eval=range(tmax)) #solve the sir ode for S,I,R given a beta and gamma
    t =  sir_sol['t']
    S = sir_sol['y'][0]
    I = sir_sol['y'][1]
    R = sir_sol['y'][2]
    return (t,S,I,R)


def plot_results(results,C):
    fig, ax = plt.subplots(1,3)
    tmax = 2*len(C)
    (t, S, I, R) = solveode((results['S'],results['I0'],results['R0']),(results['beta'],results['gamma'],tmax))

    ax[0].plot(t,S)
    ax[0].set_title('S')
    ax[1].plot(t,I)
    ax[1].set_title('I')
    ax[2].plot(t,R)
    ax[2].set_title('R')
    plt.show()

#main script

country = 'china'
C=loadData(country)
results = parmest(C)
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(results)
plot_results(results,C)



