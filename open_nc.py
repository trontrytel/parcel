from scipy.io import netcdf
import numpy as np
import time
import datetime as dt
import matplotlib.pyplot as plt



def acnvUpdate(ctY, rtY, atY, cxVs, rxVs, axVs, cdxVs, rdxVs, adxVs, j, prev, acnv):
 """ Updates the autoconversion rate matrix at each timestep
 """
 rdatY = rtY[j][:]
 acnv = np.sum(rdatY)-prev

 return(acnv)



def plotData(ctempY, rtempY, atempY, cxVals, rxVals, axVals, cdxVals, rdxVals, adxVals, i, stamp):
 """ Plots the value for aerosol, cloud, and rain at each timestep
 """
 # creates line that seperates cloud vs. rain  
 yMax = 10**-7
 lineX = [cxVals[-1] + cdxVals[-1], cxVals[-1] + cdxVals[-1]]
 lineY = [0, yMax]
 lineX2 = [axVals[-1] + adxVals[-1], axVals[-1] + adxVals[-1]]
 lineY2 = [0, yMax]

 # plots aerosol
 #plt.xlim(, rxVals[-1] + rdxVals[-1])
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

 plt.title('Cloud to Rain Autoconversion\n')
 plt.xlabel('wet radius ($m$)')
 plt.ylabel('$3^{rd}$ moment')
 plt.yscale('log')
 plt.xscale('log')
 plt.pause(0.25)
 plt.savefig('../../My_Files/figures/' + stamp + '_plot.pdf')
 plt.clf()
plt.show()



def open_nc(dct):
 """ Choose the moment to plot and the radii type; writes the data to a text file
 """
 # access to netCDF file
 f = netcdf.netcdf_file('test_acnv.nc', 'r', mmap=False)
 ts = time.time()
 stamp = dt.datetime.fromtimestamp(ts).strftime("%Y_%m_%d__%H_%M_%S")
 file = open('../../My_Files/acnv_data/' + stamp + '_acnv_data.txt', 'w')

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

 # writes the dictionary values to the file
 file.write(str(dct))

 # writes each array of variables to the file
 for key in f.variables:
  temp = f.variables[key]
  try:
   iMax, jMax = np.array(temp[:]).shape                                                                                                 
  except ValueError:
   iMax = len(temp[:])                                                                                                                  
   jMax = 0
  file.write("\n" + str(key) + ", " + str(temp) + "\n")
  file.write(str(temp[:]) + "," + "\n")

 # creates and fills the array acnv_rate for each dt
 acnv_rates = np.zeros(iMax)
 r_prev = 0

 for i in range(iMax):
  acnv_rates[i] = acnvUpdate(ctempY, rtempY, atempY, cxVals, rxVals, axVals, cdxVals, rdxVals, adxVals, i, r_prev, acnv_rates)
  plotData(ctempY, rtempY, atempY, cxVals, rxVals, axVals, cdxVals, rdxVals, adxVals, i, stamp)
  r_prev = acnv_rates[i]

 file.write("\nacnv rates: \n" + str(acnv_rates))
 file.close()
