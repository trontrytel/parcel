from ensembleAverage import ensembleAverage
from getError import getError
import matplotlib.lines as mlines
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os


def getFitAndErr(pDict, Xs, folder, varyVar, par):
    yVals = {}
    errVals = {}

    for i in range(len(pDict)):
        yVals[i]=[]
        errVals[i]=[]

    disp = ['0.01', '0.04']; n = 0
    dirPath = os.path.dirname(os.path.realpath(__file__))
    avgsDict = ensembleAverage(str(dirPath)+str(folder), str(varyVar), Xs)
    for i in range(len(avgsDict)//5):
        try: 
            Nds = avgsDict['Nd '+str(Xs[i])]
            colls = avgsDict['acnv '+str(Xs[i])]
            qRs = avgsDict['qR '+str(Xs[i])]
            rhos = avgsDict['rho '+str(Xs[i])]
            qCs = avgsDict['qC '+str(Xs[i])]
        except KeyError:
            Nds = avgsDict['Nd '+str(Xs[i])+disp[n]]
            colls = avgsDict['acnv '+str(Xs[i])+disp[n]]
            qRs = avgsDict['qR '+str(Xs[i])+disp[n]]
            rhos = avgsDict['rho '+str(Xs[i])+disp[n]]
            qCs = avgsDict['qC '+str(Xs[i])+disp[n]]
            n += 1

        fit, err = getError(Nds, qRs, qCs, rhos, colls, par)
        for j in range(len(pDict)):
            yVals[j].append(fit[j])
            errVals[j].append(err[j])
    return yVals, errVals


# input name of parameterization, returns out ordered parameters
def getParams(person):
    d = {}
    if person == 'Kessler':
        d = {0:'B', 1:'$q_{C,0}]$', 2:'a', 3:'b', 4:'c'}
    elif person == 'Beheng':
        d = {0:'C', 1:'$d_{1}$', 2:'$d_{2}$', 3:'a', 4:'b', 5:'c', 6:'e'}
    elif person == 'Khairoutdinov and Kogan':
        d = {0:'A', 1:'a', 2:'b', 3: 'd', 4:'e'}
    elif person == 'Tripoli and Cotton':
        d = {0:'D', 1:'a', 2:'b', 3:'$q_{C,0}$', 4:'c'}
    return d


def plotParamsAndErrs(pDict, xVals, valDictY, valDictErr, xlab, name, varyParam, allX, allLab):
    colorDict = {'w':'r', 'coal_kernel':'orange', 'SD':'g', 'gccn':'b'}
    
    for i in range(len(pDict)):
        errVals = valDictErr[i]
        yVals = valDictY[i]
        equConst = pDict[i]
        fitVal = pDict[i]

        # plots the fit and error for each param in the parameterization
        plt.figure(i+2)
        size = 'x-large'
        a = plt.subplot(111)
        a.errorbar(xVals[:len(yVals)], yVals, yerr=errVals, ls='none', marker='o', mfc=colorDict[varyParam], ecolor=colorDict[varyParam], mec=colorDict[varyParam], markersize=8, capsize=8, xerr=None)
        a.set_title('Best Fit for const. '+str(equConst)+' in ' + str(name), fontsize=size)
        a.set_xlabel('varied parameter', fontsize=size)
        a.set_ylabel(str(fitVal), fontsize=size)
        a.tick_params(axis='both', labelsize='large')

        # legend and formatting
        w = mlines.Line2D([],[],color=colorDict['w'], label='velocity [m/s]')
        coal_kernel = mlines.Line2D([],[],color=colorDict['coal_kernel'], label='collision kernel')
        SD = mlines.Line2D([],[],color=colorDict['SD'], label='aerosol standard deviation')
        gccn = mlines.Line2D([],[],color=colorDict['gccn'], label='giant aerosols')
        plt.legend(handles=[w, coal_kernel, SD, gccn])
        locs, labels = plt.xticks(allX, allLab)
        plt.setp(labels, rotation=30)
        a.set_xlabel('varied parameter', fontsize=size)


if __name__ == '__main__':
    """param: can only be 'Khairoutdinov and Kogan', 'Kessler', 'Beheng', or \
              'Tripoli and Cotton'
       folder : where the files to iterate over are located
       varyPar: the name (according to parcel.py model) which you are varying\
                (ex. 'w', 'coal_kernel', 'SD', 'gccn')
       vary: title for the x-axis label for the parameter the user varies
       strArr : array of strings for the titles (according to parcel.py model) \
                that user is varying (ex. strArr=['onishi_hall', 'hall_pinsky_stratocumulus',\
                'hall_pinsky_cumulonimbus', 'hall'])
    """

#*******************************************************************************************************************************************#

    varyDict = {'gccn':'giant aerosols', \
                'SD':'SD', \
                'coal_kernel':'collision kernel', \
                'w':'velocity [m/s]'}

    strDict = {'gccn':['{"kappa" : 1.28, "r_d" : [10e-7, 80e-7], "n_tot" : [8, 2]}', '{"kappa" : 1.28, "r_d" : [10e-7, 80e-7], "n_tot" : [10, 4]}', '{"kappa" : 1.28, "r_d" : [10e-7, 80e-7], "n_tot" : [15, 5]}', '{"kappa" : 1.28, "r_d" : [10e-7, 80e-7], "n_tot" : [20, 8]}', '{"kappa": 1.28, "r_d" : [10e-7, 80e-7], "n_tot" : [50, 80]}', '{"kappa" : 1.28, "r_d" : [10e-7, 80e-7], "n_tot" : [100, 60]}'], \
               'SD':['{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.026e-6], "gstdev": [1.5], "n_tot": [267e6]}}', '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.026e-6], "gstdev": [1.7], "n_tot": [267e6]}}', '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.026e-6], "gstdev": [1.9], "n_tot": [267e6]}}', '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.026e-6], "gstdev": [2.1], "n_tot": [267e6]}}'],\
               'coal_kernel':['hall_pinsky_stratocumulus', 'hall_pinsky_cumulonimbus', 'hall_pinsky_stratocumulus', 'hall', 'onishi_hall', 'onishi_hall'],\
               'w':[0.75, 0.25, 0.5, 0.75, 1.0]}

    labDict = {'gccn':['[100, 60]', '[15, 5]', '[10, 4]', '[20, 8]', '[50, 80]', '[8, 2]'],\
               'SD':['1.5', '1.7', '1.9', '2.1'],\
               'coal_kernel':['HPS', 'HPS', 'H', 'OH1', 'OH4'],\
               'w':['0.25', '0.5', '0.75', '1.0']}

    folderDict = {'gccn':'/vary_gccns/',\
                  'SD':'/vary_SD/',\
                  'coal_kernel':'/vary_kerns_final/',\
                  'w':'/vary_vel3/'}


    # input the parameterization and the array of variables to iterate over
    param = 'Kessler'
    varyPars = ['coal_kernel']

    allXs = []; allLabs = []
    n = 1
    for i in range(len(varyPars)):
        folder = folderDict[varyPars[i]]
        vary = varyDict[varyPars[i]]
        strArr = strDict[varyPars[i]]
        labs = labDict[varyPars[i]]
        xArr = np.linspace(n, len(labs)+n-1, len(labs))
        n += len(labs)
        #*******************************************************************************************************************************************#

        paramDict = getParams(param)
        fits, errs = getFitAndErr(paramDict, strArr, folder, varyPars[i], param)

        plotParamsAndErrs(paramDict, xArr, fits, errs, vary, param, varyPars[i], xArr, labs)

    plt.show()
