from scipy.io import netcdf
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os


# finds the index at which the temperature array falls below 273.15K (0C)
def getTempEndInd(arr):
    for i in range(len(arr)):
        if arr[i] <= 273.15:
            return i
    return len(arr)


# finds the index at which the cloud has decreased to 0.25 times its max number of droplets
def getCloudEndInd(arr):
    ind = 0; end = 0
    ind = np.argmax(arr)
    for j in range(ind+1, len(arr)):
        if arr[j]<=0.25*arr[ind]:
            return j
    return len(arr)


# returns raw data from the .nc file given as an input
def getData(fName):
    f = netcdf.netcdf_file(fName, 'r', mmap=False)
    dataset = Dataset(fName)
    vel = 0; coal_kernel = ''; z_max = 0; aerosol = ''; diss = 0; giantAero = ''

    # handling the attributes in the file
    for attr in dataset.ncattrs():
        if attr ==  'w':
            vel = float(getattr(dataset, attr))
        elif attr == 'coal_kernel':
            coal_kern = str(getattr(dataset, attr))
        elif attr == 'z_max':
            z_max = float(getattr(dataset, attr))
        elif attr == 'aerosol':
            aero = str(getattr(dataset, attr))
        elif attr == 'coal_dissipation_rate':
            diss = float(getattr(dataset, attr))
        elif attr == 'gccn':
            giantAero = str(getattr(dataset, attr))

    # handling the data arrays in the file
    N_d = f.variables['cradii_m0'][:]
    numR = f.variables['rradii_m0'][:]
    q_C = f.variables['cradii_m3'][:]*4./3.*np.pi*1000
    q_R = f.variables['rradii_m3'][:]*4./3.*np.pi*1000
    rho_d = f.variables['rhod'][:]
    ACNV = f.variables['acnv'][:]
    T = f.variables['T'][:]
    t = f.variables['t'][:]
    condr = f.variables['cond_growth_rain'][:]

    return N_d, numR, q_C, q_R, ACNV, vel, coal_kern, aero, z_max, T, diss, rho_d, giantAero, t, condr


# returns the length of the acnv (all collisions) array
def arrLens(fName):
    f = netcdf.netcdf_file(fName, 'r', mmap=False)
    T = f.variables['T'][:]
    return len(T)


# returns the data from all files concatenated and sorted by temperature
def getSortedArrs(direc, fileArr):
    Nds, Nrs, qCs, qRs, acnvs, ws, coal_kernels, aerosols, z_maxs, Ts, disps, rhos, gccn, time, cond = getData(direc+fileArr[0])
    indT = getTempEndInd(Ts)
    indNd = getCloudEndInd(Nds)

    end = len(rhos)
    if indT < indNd:
        end = indT
    else:
        end = indNd

    # ends the arrays where they are no longer valid for the first file
    Nds = Nds[0:end]
    acnvs = acnvs[0:end]
    qCs = qCs[0:end]
    qRs = qRs[0:end]
    rhos = rhos[0:end]
    Ts = Ts[0:end]

    endT = int(1e10); endNd = int(1e10) # starts at very large so that index compared later will be smaller
    for j in range(1, len(fileArr)):
        Nd, Nr, qC, qR, acnv, w, coal_kernel, aerosol, z_max, T, disp, rho, gccn, time, cond = getData(direc+fileArr[j])
        acnv[np.isnan(acnv)] = 0. # since 'nan' if no collisions happen; change to 0 
        indT = getTempEndInd(T)
        indNd = getCloudEndInd(Nd)
        
        if indT < indNd:
            end = indT
        else:
            end = indNd
            
        # ends the arrays where they are no longer valid
        Nd = Nd[0:end]
        acnv = acnv[0:end]
        qC = qC[0:end]
        qR = qR[0:end]
        rho = rho[0:end]
        T = T[0:end]

        # concatenates all of the respective arrays into one
        Nds = np.concatenate([Nds, Nd], axis=0)
        acnvs = np.concatenate([acnvs, acnv], axis=0)
        qCs = np.concatenate([qCs, qC], axis=0)
        qRs = np.concatenate([qRs, qR], axis=0)
        rhos = np.concatenate([rhos, rho], axis=0)
        Ts = np.concatenate([Ts, T], axis=0)

    # orders the arrays all relative to T (corresponds to height)
    order = np.argsort(Ts, kind='mergesort', axis=0)
    Nds = Nds[order]
    acnvs = acnvs[order]
    qRs = qRs[order]
    rhos = rhos[order]
    qCs = qCs[order]
    return Nds, acnvs, qCs, qRs, rhos


# plots collision versus condensation; points are colored by height
def ScatterACNV(fName, j):
    size = 'x-large'
    plt.figure(j)
    fig = plt.subplot(111)
    fig.set_xlabel('rain condensation $[kg  s^{-1}  m^{-3}]$', fontsize=size)
    fig.set_ylabel('collision rate $[kg  s^{-1}  m^{-3}]$', fontsize=size)
    fig.set_xscale('log')
    fig.set_yscale('log')
    fig.set_xlim(1e-10, 1e-5)
    fig.set_ylim(1e-9, 1e-4)
    fig.tick_params(axis='both', labelsize='large')

    Nd, Nr, qC, qR, acnv, w, coal_kernel, aerosol, z_max, T, disp, rho, gccn, time, cond_r  = getData(fName)

    indC = getCloudEndInd(Nd)
    indT = getTempEndInd(T)

    end = 0
    if indC > indT:
        end = indT
    else:
        end = indC

    h = time[:]*w; ratios=h[:]/h[-1]
    fig.set_title('collision rate vs. rain condensation '+str(coal_kernel)+' at '+str(w)+' m/s', fontsize=size)
    ratios = ratios[acnv>1e-9]
    cond_r = cond_r[acnv>1e-9]
    acnv = acnv[acnv>1e-9]
    plt.scatter(cond_r[:end], acnv[:end], c=ratios[:end])
    plt.colorbar()


def ensembleAverage(directory, var, vals):
    """ directory: the directory containing the desired files
        var: the variable (coal_kernel, w, aerosol) being looked at
        vals: the associated options for the given variable (ex. var='w', vals=[1.0, 1.5, 2.0]

        Uses the var to return a dictionary for values of yVal and errVal for each ensemble average
        for each value in the vals array. Vals can then be compared to one another.
    """
    files = os.listdir(directory)
    tol = 0.001
    
    # creates a dictionary that associates vals with the file names
    nameDict = {}
    for i in range(len(vals)):
        if vals[i] == 'onishi_hall':
            nameDict[vals[i]+'0.01'] = []
            nameDict[vals[i]+'0.04'] = []
        else:
            nameDict[vals[i]] = []
    
    # groups together the files that all have the same val[i]
    #     (ex. if var='w', files are grouped by each velocity)
    for filename in files:
        dataDict = {}
        if filename.endswith('.nc'):
            Nd, Nr, qC, qR, acnv, w, coal_kernel, SD, z_max, T, disp, rho, gccn, t, cond = getData(directory+filename)

            dataDict['w']=w
            dataDict['coal_kernel']=coal_kernel
            dataDict['SD']=SD
            dataDict['gccn']=gccn

            try:
                for i in range(len(vals)):
                    vals[i] = float(vals[i])
                    if abs(dataDict[var]-vals[i]) <= tol:
                        nameDict[vals[i]].append(filename) # appends the filename to the dictionary with its val

            except ValueError:
                if var == 'coal_kernel':
                    for i in range(len(vals)):
                        if dataDict[var]==vals[i]: # references the specific 'var' (input) in dataDict
                            if vals[i] == 'onishi_hall':

                                if abs(disp-0.04) <= tol:
                                    num = '0.04'
                                    nameDict[vals[i]+num].append(filename)
                                else:
                                    num = '0.01'
                                    nameDict[vals[i]+num].append(filename)
                            else:
                                nameDict[vals[i]].append(filename) # appends the filename to the dictionary with its val

                else:
                    for i in range(len(vals)):
                        if dataDict[var]==vals[i]: # references the specific 'var' (input) in dataDict 
                            nameDict[vals[i]].append(filename) # appends the filename to the dictionary with its val

    # creates dictionary of sorted arrays corresponding to their 'val'
    outDict = {}
    n = 0
    disip = ['0.01', '0.04']
    for i in range(len(vals)):
        if vals[i]=='onishi_hall':
            files = nameDict[vals[i]+disip[n]]
            Nd, acnv, qC, qR, rho = getSortedArrs(directory, files)
            
            outDict['Nd '+str(vals[i])+disip[n]] = Nd
            outDict['acnv '+str(vals[i])+disip[n]] = acnv
            outDict['qC '+str(vals[i])+disip[n]] = qC
            outDict['qR '+str(vals[i])+disip[n]] = qR
            outDict['rho '+str(vals[i])+disip[n]] = rho
            n += 1

        else:
            files = nameDict[vals[i]]
            Nd, acnv, qC, qR, rho = getSortedArrs(directory, files)

            outDict['Nd '+str(vals[i])] = Nd
            outDict['acnv '+str(vals[i])] = acnv
            outDict['qC '+str(vals[i])] = qC
            outDict['qR '+str(vals[i])] = qR
            outDict['rho '+str(vals[i])] = rho

    return outDict
