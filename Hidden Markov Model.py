import calendar
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib import style
import os
import pandas as pd

from scipy.stats import norm
import numpy as np
from numpy.linalg import inv
from scipy import linalg

style.use('ggplot')

iStates = 2
iMaxiter = 400
nTol = 10e-7
data=pd.read_csv('C:\AMS 518\ETFs.csv')
vnData=data['closeunadj'].tolist()


def xMarkovChainEquilibrium(mnTransitionMatrix):
    a = np.dot([linalg.inv(mnTransitionMatrix.transpose() - np.eye(len(mnTransitionMatrix)) + 1)],(constArray(mnTransitionMatrix)))
    # print(a)
    return a

def xHMMUGRandomInitial(vnData, iStates):
    global mnA
    global vnIota
    global vnM
    global vnS
    mnA = (np.random.uniform(0.01, 0.99, size=(iStates, iStates)) / sum(np.random.uniform(0.01, 0.99, size=(iStates, iStates))))
    vnIota = xMarkovChainEquilibrium(mnA).flatten()
    vnM = ((np.random.uniform(0.01, 0.99, size=iStates)) / sum(np.random.uniform(0.01, 0.99, size=iStates))) * np.mean(
        vnData)
    vnS = np.sqrt(
        ((np.random.uniform(0.01, 0.99, size=iStates)) / sum(np.random.uniform(0.01, 0.99, size=iStates))) * np.var(
            vnData))
    print('These are the values of A, Iota, M, S after random initialize')
    print('@@@@@@@Q@@@@@@@@@@')
    print('mnA')
    print(mnA)
    v1 = mnA
    print('vnIota')
    print(vnIota)
    v2 = vnIota
    print('vnM')
    print(vnM)
    v3 = vnM
    print('vnS')
    print(vnS)
    v4 = vnS
    return mnA, vnIota, vnM, vnS


def xHHMUG(vnData, vnInitialIota, mnInitialLA, vnInitialMean, vnInitialSdev):
    print('Here are values of iT, iN, vnIota, mnA, vnMean, vnSdev at start of function')
    iT = len(vnData)
    print('iT is:')
    print(iT)
    iN = len(mnInitialLA)
    print('iN is:')
    print(iN)
    vnIota = vnInitialIota
    mnA = mnInitialLA
    vnMean = vnInitialMean
    vnSdev = vnInitialSdev
    mnB = []
    for x in range(iT):
        for y in range(iN):
            mnB.append(norm.pdf(vnData[x], loc=vnMean[y], scale=vnSdev[y]))
    print('mnB is:')
    print(mnB)
    # mnAlpha = np.empty((iT, iN))
    mnAlpha = np.zeros((iT, iN))
    print('mnAlpha is pt 1:')
    print(mnAlpha)
    vnNu = np.zeros(iT)
    print('vnNu is pt 1:')
    print(vnNu)
    # vnNu = np.empty(iT)
    mnAlpha[0] = np.multiply(vnIota, mnB[0])
    print('mnAlpha after zero pt 1')
    print(mnAlpha)
    vnNu[0] = np.divide(1, np.sum(mnAlpha[0]))
    print('vnNU after zero pt 1')
    print(vnNu)
    mnAlpha[0] = np.multiply(mnAlpha[0], vnNu[0])
    print('mnAlpha[0] after zero pt 2')
    print(mnAlpha[0])
    for x in range(1, iT):
        # first mnAlpha is changed to dot from multiply
        mnAlpha[x] = np.multiply(np.dot(mnAlpha[x - 1], mnA), mnB[x])
        vnNu[x] = np.divide(1, np.sum(mnAlpha[x]))
        mnAlpha[x] = np.multiply(mnAlpha[x], vnNu[x])
    print('values after for loop')
    print('alpha')
    print(mnAlpha)
    print('Nu')
    print(vnNu)
    print('alpha again')
    print(mnAlpha)
    # vnLogLikelihood = (-np.sum(np.log(vnNu)))
    vnLogLikelihood = []
    print('vnlog')
    print(vnLogLikelihood)
    ###################### E step ###############################

    print('Nu E step')
    print(vnNu)
    vnLogLikelihood.append(-np.sum(np.log(vnNu)))
    i = 0

    while i < iMaxiter:
        mnBeta = np.zeros((iT, iN))
        # vnNu = np.empty(iT)
        mnBeta[iT-1] = vnNu[iT-1]
        k = iT - 2
        while k >= 0:
            mnBeta[k] = np.dot(mnA, (np.multiply(mnBeta[k+1], mnB[k+1])*vnNu[k]))
            k = k - 1



        # for x in range(iT - 1, 0):
        #     mnBeta[x] = np.multiply(np.dot(mnA, np.multiply(mnBeta[x + 1], mnB[x + 1])), vnNu[x])
        print('mnBeta')
        print(mnBeta)
        mnGamma = np.multiply(mnAlpha, mnBeta)
        # mnGamma = np.divide(np.multiply(mnAlpha, mnBeta), np.sum(np.multiply(mnAlpha, mnBeta)))
        for i in range(len(mnGamma)):
            mnGamma[i] = mnGamma[i]/np.sum(mnGamma[i])
        print('gamma')
        # print('gama')
        print(mnGamma)
        # mnXi = []
        mnXi = np.empty((iN, iN))
        for x in range(iT - 1, 0):
            p = mnAlpha[x]
            q = np.multiply(mnBeta[x + 1], mnB[x + 1])
            a = np.multiply(mnA, np.kron(p, q))
            k = np.divide(a, np.sum(a.flatten()))
            np.append(mnXi, k)
            # mnXi.append(np.divide(a, np.sum(a.flatten())))
        print('mnXi')
        print(mnXi)
        ########################## M-step ################################
        mnA = (mnXi / np.sum(mnXi))
        print('mnA')
        print(mnA)

        vnIota = mnGamma[-1]
        mnWeights = np.divide(mnGamma.transpose(), np.sum(mnGamma.transpose()))
        print('weights')
        print(mnWeights)
        # change below line to iN to see diff results - debugging notes ignore
        vnMean = np.dot(mnWeights, vnData)
        print('mean')
        print(vnMean)

        vnSdev= [0]*iN

        for x in range(iN):
            ex = np.array((vnData-vnMean[x])**2)
            vnSdev[x] = np.sqrt(np.sum(np.multiply(ex, mnWeights)))
        print('sdList')
        print(vnSdev)

        mnB = []
        for x in range(iT):
            for y in range(iN):
                mnB.append(norm.pdf(vnData[x], loc=vnMean[y], scale=vnSdev[y]))

        print('this is initial mnB')
        print(mnB)
        # changed np.empty to np.zeros
        mnAlpha = np.zeros((iT, iN))

        vnNu = np.zeros(iT)

        mnAlpha[0] = np.multiply(vnIota, mnB[0])

        vnNu[0] = np.divide(1, np.sum(mnAlpha[0]))

        mnAlpha[0] = np.multiply(mnAlpha[0], vnNu[0])
        print('alpha M step')

        print(mnAlpha)
        for x in range(1, iT):
            mnAlpha[x] = np.multiply(np.dot(mnAlpha[x - 1], mnA), mnB[x])

            vnNu[x] = np.divide(1, np.sum(mnAlpha[x]))
            # mnAlpha[x] *= vnNu[x]
            mnAlpha[x] = np.multiply(mnAlpha[x], vnNu[x])

        # vnLogLikelihood = []
        print('Nu M step')
        print(vnNu)
        print('len of Nu')
        print(len(vnNu))
        # for i in range(len(vnNu)):
        #     logs = np.log(vnNu)
        #     vnLogLikelihood[i] = -1*(np.sum(logs))
            # vnLogLikelihood[i] = (-np.sum(np.log(vnNu)))
        # vnLogLikelihood.append((-np.sum(np.log(vnNu))))
        # np.append([vnLogLikelihood], [((-np.sum(np.log(vnNu))))])
        # print('logs value')
        # print(logs)


        # vnLogLikelihood = (-np.sum(np.log(vnNu)))
        print('loglike2')
        print(vnLogLikelihood)
        vnLogLikelihood.append(-np.sum(np.log(vnNu)))
        i = i + 1
        if (vnLogLikelihood[-1] / (vnLogLikelihood[-2] - 1)) <= nTol:
            break
    print('iN')
    print(iN)
    v = (iN*(iN+2)-1)

    print('nBIC')

    nBIC = -2 * vnLogLikelihood[-1] + v * np.log(iT)

    print('this is the results')
    print(nBIC)
    print('vnIota ')
    print(vnIota)
    print(' mnA ')
    print(mnA)
    print('vnMean')
    print(vnMean)
    print('vnSdev')
    print(vnSdev)
    print('mnGamma')
    print(mnGamma)
    print('BIC')
    print(nBIC)
    print('logLike')
    print(vnLogLikelihood)

    return vnIota, mnA, vnMean, vnSdev, mnGamma, nBIC, vnLogLikelihood


