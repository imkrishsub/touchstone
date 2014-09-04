#!/usr/bin/python
import numpy
#import math
import matplotlib.pyplot as plt
import os.path

from SchoofSemiAnalytic import schoofSemiAnalyticA, schoofSemiAnalyticB

Nx = 1025
As = [30,
      28.231,
      26.566,
      25,
      23.210,
      21.544,
      20,
      18.171,
      16.510,
      15,
      13.103,
      11.447,
      10,
      7.9370,
      6.2996,
      5,
      3.9685,
      3.1498,
      2.5,
      2.3208,
      2.1544,
      2,
      1.8171,
      1.6510,
      1.5,
      1.3104,
      1.1447,
      1,
      0.79279,
      0.62996,
      0.5,
      0.39685,
      0.31498,
      0.25]
As = numpy.array(As)*1e-26 #Pa^(-3) s^(-1)



folders = []
folders.append('./dense_Ab/retreated_Nx_%i' % (Nx))
folders.append('./dense_Ab/advanced_Nx_%i' % (Nx))
folders.append('./dense_Ab/unstable_Nx_%i' % (Nx))

colors = ['r.', 'm.', 'b.','c.','g.', 'y.']

ps = numpy.linspace(0.0,1.0,5)

xgMin = 1.5
xgMax = 0
for AIndex in range(len(As)):
  pIndex = 0
  p = ps[pIndex]    
  A = As[AIndex]
  for regionIndex in range(3):
    fileName = '%s/p_%.4f_A_%.4e.da'%(folders[regionIndex],p,A)
    if(not os.path.exists(fileName)):
      continue
    filePointer = open(fileName,'rb')
    filePointer.seek(4*8*Nx,0) # skip H, u and x
    xg = numpy.fromfile(filePointer, dtype=float, count=1)
    xgMax = max(xgMax,xg)
    xgMin = min(xgMin,xg)

xgs = numpy.linspace(xgMin,xgMax,1025)
AschoofA = schoofSemiAnalyticA(xgs)
AschoofB = schoofSemiAnalyticB(xgs)
plt.semilogx(1.0/AschoofA,xgs,'k')
plt.semilogx(1.0/AschoofB,xgs,'r')

plt.figure(1)
for AIndex in range(len(As)):
  for pIndex in range(len(ps)):
    p = ps[pIndex]    
    A = As[AIndex]
    for regionIndex in range(3):
      fileName = '%s/p_%.4f_A_%.4e.da'%(folders[regionIndex],p,A)
      if(not os.path.exists(fileName)):
        continue
      filePointer = open(fileName,'rb')
      filePointer.seek(4*8*Nx,0) # skip H, u and x
      xg = numpy.fromfile(filePointer, dtype=float, count=1)
      plt.semilogx(1.0/A,xg,colors[pIndex])

labels = []
labels.append('Schoof 2007 model A')
labels.append('Schoof 2007 model B')
for pIndex in range(len(ps)):
  labels.append('p = %.2f'%ps[pIndex])

plt.legend(labels, loc='lower right')

plt.show()

# vim: ts=2 noet
