import numpy as np
import copy
import operator
import scipy.optimize as sciopt
from scipy.io import netcdf
import matplotlib.pyplot as plt

np.set_printoptions(precision=40)

# Khairoutdinov & Kogan (2000)        #*****Doesn't include rho!!!*****
def KhairKogFunc(x, A, a, b, d, e):
    Nd  = x[0];     qC   = x[1];     qR = x[2]
    rho = x[3];     coll = x[4]
    outArr = A*qC[:]**a*Nd[:]**(b)+d*(qC[:]*qR[:])**e
    return outArr

#  Kessler (1969)
def KesslerFunc(x, B, qC0, a, b, c):
    Nd  = x[0];    qC   = x[1];   qR = x[2]
    rho = x[3];    coll = x[4]
    maxArr = qC[:]-qC0
    maxArr[maxArr<0]=0.
    outArr = B*maxArr[:]+a*qC[:]*(qR[:])**(b)*(Nd[:])**(c)
    return outArr

# Beheng (1994)
def BehengFunc(x, C, D1, D2, a, b, c, d):
    Nd  = x[0];     qC   = x[1];     qR  = x[2]
    rho = x[3];     coll = x[4]
    D = copy.copy(Nd[i])
    D[D<200*100**3]=D1
    D[D>=200*100**3]=D2
    outArr = C*D[:]**(a)*qC[:]**(b)*Nd[:]**(c)+d*qR[:]*qC[:]
    return outArr

# Tripoli & Cotton (1980)
def TripCottFunc(x, D, a, b, qC0, c):
    Nd  = x[0];     qC   = x[1];      qR  = x[2]
    rho = x[3];     coll = x[4]
    val = qC[i]-4./3.*np.pi*1000*qC0*Nd[:]
    ind = val[val==0]
    val[val<0]=0.
    val[val>0]=1.
    val[ind]=0.5
    outArr = D*qC[:]**(a)*Nd[:]**(b)*val[:]+c*qC[:]*qR[:]
    return outArr

# averages array into 'num' number of sections
def avgArr(arr, num):
    outArr = np.zeros(len(arr)//num)
    for i in range(len(arr)//num):
        if len(arr)-i*num >= num:
            outArr[i] =sum(arr[i*num:(i+1)*num])/num
        else:
            outArr[i] =sum(arr[i*num:len(arr)])/(len(arr)-i*num)
    return outArr

# returns five averaged arrays in order
def avgVals(arr1, arr2, arr3, arr4, arr5):
    n = 30
    arr1 = avgArr(arr1, n)
    arr2 = avgArr(arr2, n)
    arr3 = avgArr(arr3, n)
    arr4 = avgArr(arr4, n)
    arr5 = avgArr(arr5, n)
    return arr1, arr2, arr3, arr4, arr5

# plots the values for collision vs qC and the best fit
def plotFit(actIn, actColl, fitParams, name, func):
    c1 = fitParams[0]; c2 = fitParams[1]; c3 = fitParams[2]; c4 = fitParams[3]
    c5 = fitParams[4]
    qC = actIn[1]*1000
    coll = actIn[4]

    size='x-large'
    fitY = func(actIn, c1, c2, c3, c4, c5)
    a = plt.subplot(111)
    a.plot(qC[:]*1000, coll[:]*10e6, c='b', marker='o', linestyle='None')
    a.plot(x, fitY[:]*10e6, linestyle='-.', color='k')
    a.set_xlabel('$q_{C}$ [g/kg]', fontsize=size)
    a.set_ylabel('Collision Rate ($10^{-6}$) [$kg$ $s^{-1}$ $m^{-3}$]', fontsize=size)
    a.set_title(name+' Best Fit', fontsize=size)
    a.tick_params(axis='both', labelsize='large')




def getError(Nd, qR, qC, rho, collRate, name):
    """Nd: numpy array for the number of cloud droplets at each time step
       collRate: numpy array for the collisions (acnv+acc) rate at each time step
       qR: numpy array for the rain water mixing ratio [kg rain water/kg dry air]
       at each time step
       qC: numpy array for the cloud water mixing ratio [kg cloud water/kg dry
       air] at each time step
       rho: numpy array for the dry air density at each time step
       name: string for the name of the parameterization desired ('Kessler', 'Khairourdinov and Kogan', etc.)
    """

    func     = {'Khairoutdinov and Kogan':KhairKogFunc,\
                'Tripoli and Cotton':TripCottFunc,\
                'Kessler':KesslerFunc,\
                'Beheng':BehengFunc}

    bnd      = {'Khairoutdinov and Kogan':((0, 1, -np.inf, 0, 0),\
                (np.inf, np.inf, 0, np.inf, np.inf)),\
                'Tripoli and Cotton':((0, 1, -1, 1e-6, 0), (np.inf, np.inf, 0, 20e-6, np.inf)),\
                'Kessler':((0, 0, 0, 0, 0), (1, 1, 1, 1, 1)),\
                'Beheng':((0, 0, 0, -np.inf, 0, -np.inf, 0), (np.inf, 15,\
                np.inf, 0, np.inf, 0, np.inf))}

    init     = {'Khairoutdinov and Kogan':[7.42e13, 2.47, -1.79, 67., 1.15],\
                'Tripoli and Cotton':[3268, 7./3., -1./3., 7e-6, 4.7],\
                'Kessler':[1e-3, 5e-4, 0.29, 7./8., 1./8.],\
                'Beheng':[3.0e34, 9.9, 3.9, -1.7, 4.7, -3.3, 6.0]}
    
    # force array sizes to be the same
    Nd = np.squeeze(Nd)
    qR = np.squeeze(qR)
    qC = np.squeeze(qC)
    rho = np.squeeze(rho)
    collRate = np.squeeze(collRate)

    # parameterizations don't work for Nd=0
    qC = qC[Nd>0]
    qR = qR[Nd>0]
    rho = rho[Nd>0]
    collRate = collRate[Nd>0]
    Nd = Nd[Nd>0]

    # *************************************************************************************** #
    fits, cov = sciopt.curve_fit(func[name], [Nd, qC, qR, rho, collRate], collRate, \
               p0=init[name], bounds=bnd[name], sigma=None, maxfev=5000)
    
    error = np.zeros(len(cov))
    for i in range(len(fits)):
        try:
            error[i] = np.sqrt(np.diag(cov))[i]
        except:
            error[i] = 0
    # *************************************************************************************** # 
    
    #NdAvg, qCAvg, qRAvg, rhoAvg, collAvg = avgVals(Nd, qC, qR, rho, collRate)
    #plotFit([NdAvg, qCAvg, qRAvg, rhoAvg, collAvg], collAvg, fits, name, func[name])

    return fits, error
