from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

np.set_printoptions(precision=25)
np.set_printoptions(threshold=np.inf)


def acnvUpdate(ctY, rtY, atY, cxVs, rxVs, axVs, cdxVs, rdxVs, adxVs, j, prevR, prevC):
 """ Updates the autoconversion rate matrix at each timestep
 """
 delRain = rtY[j]-prevR
 delCloud = ctY[j]-prevC

 return(delRain, delCloud)



def conserveMass(v, aero, cl, rn):
 ind = len(v)-1
 init = v[0] + 4./3.*np.pi*1000*aero[0] + 4./3.*np.pi*1000*cl[0] + 4./3.*np.pi*1000*rn[0]
 final =v[ind] + 4./3.*np.pi*1000*aero[ind] + 4./3.*np.pi*1000*cl[ind] + 4./3.*np.pi*1000*rn[ind]
 
 return(final-init)


def plotData(ctempY, rtempY, atempY, cxVals, rxVals, axVals, cdxVals, rdxVals, adxVals, i, outfile, z, r_v, t, acnv):
 """ Plots the value for aerosol, cloud, and rain at each timestep
 """
 # creates line that seperates cloud vs. rain  
 yMax = 10**-7
 lineX = [cxVals[-1] + cdxVals[-1], cxVals[-1] + cdxVals[-1]]
 lineY = [0, yMax]
 lineX2 = [axVals[-1] + adxVals[-1], axVals[-1] + adxVals[-1]]
 lineY2 = [0, yMax]

 plt.figure(1)
 plt.subplot(331)
 # plots aerosol
 plt.ylim(10e-25, yMax+yMax*0.1)
 try:
  plt.bar(axVals, atempY[i], adxVals, align='edge', color='yellow')
 except ValueError:
  print("atempY: " + str(atempY[i]))
  pass
 # plots cloud
 try:
  plt.bar(cxVals, ctempY[i], cdxVals, align='edge', color='green')
 except ValueError:
  print("ctempY: " + str(ctempY[i]))
  pass
 # plots rain
 try:
  plt.bar(rxVals, rtempY[i], rdxVals, align='edge', color='blue')
 except ValueError:
  print("rtempY: " + str(rtempY[i]))
  pass

 # plots the cloud/rain seperation line                                                                                                     
 plt.plot(lineX, lineY, color='black')
 plt.plot(lineX2, lineY2, color='black')

 plt.title('Cloud to Rain Autoconversion')
 plt.xlabel('wet radius ($m$)')
 plt.ylabel('$3^{rd}$ moment')
 plt.yscale('log')
 plt.xscale('log')
 
 plt.subplot(332)
 plt.ylim(0, z[len(z)-1])
 plt.plot(r_v[0:i], z[0:i], color="blue")
 plt.title("Height versus Water Vapor")
 plt.xlabel("water vapor/dry air [$kg/kg$]")
 plt.ylabel("z [$m$]")
 
 # density of water is 1000 kg/m^3
 plt.subplot(334)
 rc = 4./3.*np.pi*1000*ctempY[:]
 plt.ylim(0, z[len(z)-1])
 plt.plot(rc[0:i], z[0:i], color="blue")
 plt.title("Height versus Cloud Water")
 plt.xlabel("cloud water/dry air [$kg/kg$]")
 plt.ylabel("z [$m$]")

 plt.subplot(335)
 rr = 4./3.*np.pi*1000*rtempY[:]
 plt.ylim(0, z[len(z)-1])
 plt.plot(rr[0:i], z[0:i], color="blue")
 plt.title("Height versus Rain Water")
 plt.xlabel("rain water/dry air [$kg/kg$]")
 plt.ylabel("z [$m$]")

 plt.subplot(336)
 ra = 4./3.*np.pi*1000*atempY[:]
 plt.ylim(0, z[len(z)-1])
 plt.plot(ra[0:i], z[0:i], color="blue")
 plt.title("Height versus Aerosol Water")
 plt.xlabel("aerosol water/dry air [$kg/kg$]")
 plt.ylabel("z [$m$]")

 plt.subplot(337)
 plt.plot(rr[0:i], acnv[0:i], color="blue")
 plt.title("Autoconversion Rate vs. Rain Water Mass")
 plt.xlabel("rain water/dry air [$kg/kg$]")
 plt.ylabel("autoconversion rate ($dq_{l}/dt$)")

 plt.subplot(338)
 plt.plot(rc[0:i], acnv[0:i], color="blue")
 plt.title("Autoconversion Rate vs. Cloud Water Mass")
 plt.xlabel("cloud water/dry air [$kg/kg$]")
 plt.ylabel("autoconversion rate ($dq_{l}/dt$)")

 plt.subplot(339)
 plt.plot(ra[0:i], acnv[0:i], color="blue")
 plt.title("Autoconversion Rate vs. Aerosol Water Mass")
 plt.xlabel("aerosol water/dry air [$kg/kg$]")
 plt.ylabel("autoconversion rate ($dq_{l}/dt$)")

 plt.subplots_adjust(left=None, bottom=0.08, right=None, top=0.96, wspace=0.3, hspace=0.3)
 plt.savefig(outfile + '_plot.pdf')


def open_nc(dct):
 """ Choose the moment to plot and the radii type; writes the data to a text file
 """
 # access to netCDF file
 outfile = dct['outfile'][0:-3]
 f = netcdf.netcdf_file(outfile + '.nc', 'r', mmap=False)

 # cloud-related constants
 cxKey = 'cradii_r_wet'
 cdxKey = 'cradii_dr_wet'
 cyKey = 'cradii_m3'
 cxVals = f.variables[cxKey][:]
 cdxVals = f.variables[cdxKey][:]
 ctempY = f.variables[cyKey]

 # rain-related constants
 rxKey = 'rradii_r_wet'
 rdxKey = 'rradii_dr_wet'
 ryKey = 'rradii_m3'
 rxVals = f.variables[rxKey][:]
 rdxVals = f.variables[rdxKey][:]
 rtempY = f.variables[ryKey]

 # aerosol-related constants
 axKey = 'aradii_r_wet'
 adxKey = 'aradii_dr_wet'
 ayKey = 'aradii_m3'
 axVals = f.variables[axKey][:]
 adxVals = f.variables[adxKey][:]
 atempY = f.variables[ayKey]

 # r_v- and z-related constants
 r_vKey = 'r_v'
 zKey = 'z'
 tKey = 't'
 r_vVals = f.variables[r_vKey][:]
 zVals = f.variables[zKey][:]
 tVals = f.variables[tKey][:]

 # writes each array of variables to the file
 for key in f.variables:
  temp = f.variables[key]
  try:
   iMax, jMax = np.array(temp[:]).shape                                                                                                 
  except ValueError:
   iMax = len(temp[:])                                                                                                                  
   jMax = 0

 # creates and fills the arrays acnv_rates and cloud_rates for each dt
 acnv_rates = np.zeros(iMax, dtype=float)
 cloud_rates = np.zeros(iMax, dtype=float)
 cArr = np.zeros(iMax, dtype=float)
 rArr = np.zeros(iMax, dtype=float)
 r_prev = 0; c_prev = 0

 for i in range(iMax):
  acnv_rates[i], cloud_rates[i] = acnvUpdate(ctempY, rtempY, atempY, cxVals, rxVals, axVals, cdxVals, rdxVals, adxVals, i, r_prev, c_prev)
  r_prev = acnv_rates[i]
  c_prev = cloud_rates[i]
 plotData(ctempY, rtempY, atempY, cxVals, rxVals, axVals, cdxVals, rdxVals, adxVals, iMax-1, outfile, zVals, r_vVals, tVals, acnv_rates)
