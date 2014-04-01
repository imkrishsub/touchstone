#!/usr/bin/python
#import numpy
from Solver import Solver
from optparse import OptionParser

#from pprint import pprint

import os.path
#import numpy.linalg

parser = OptionParser()

parser.add_option("--p", type="float", default=0.0, dest="p")
parser.add_option("--A", type="float", default=1e-24, dest="A")
parser.add_option("--C", type="float", default=7.624e6, dest="C")
parser.add_option("--W", type="float", default=10, dest="W")
parser.add_option("--rho_i", type="float", default=900.0, dest="rho_i") # kg/m^3 ice density
parser.add_option("--rho_w", type="float", default=1000.0, dest="rho_w") # kg/m^3 water density
parser.add_option("--a", type="float", default=1.0, dest="a")
parser.add_option("--g", type="float", default=9.8, dest="g")
parser.add_option("--n", type="float", default=3.0, dest="n")
parser.add_option("--linearSlope", type="float", default=778.5, dest="linearSlope") #drop in m per 750 km, as in Schoof 2007
parser.add_option("--lambda_0", type="float", default=2, dest="lambda_0")
parser.add_option("--m_0", type="float", default=2.0, dest="m_0")
parser.add_option("--Ab", type="float", default=3.1688e-24, dest="Ab") #Pa^-3 s^-1
parser.add_option("--poly", action="store_true", dest="poly")
parser.add_option("--useGLP", action="store_true", default=False, dest="useGLP")

parser.add_option("--inFile", type="string", default="none", dest="inFile")
parser.add_option("--folder", type="string", default="results", dest="folder")
parser.add_option("--outFile", type="string", default="results.pyda", dest="outFile")
parser.add_option("--filePointer", type="string", default="default.pointer", dest="filePointer")

parser.add_option("--Nx", type="int", default=1321, dest="Nx")
parser.add_option("--xc", type="float", default=2.112, dest="xc")
parser.add_option("--deltaX", type="float", default=1.6e-3, dest="deltaX") # m
parser.add_option("--maxSteps", type="int", default=20000, dest="maxSteps")
parser.add_option("--maxInnerSteps", type="int", default=200, dest="maxInnerSteps")
parser.add_option("--maxInitSteps", type="int", default=200, dest="maxInitSteps")
parser.add_option("--stepsPerWrite", type="int", default=10, dest="stepsPerWrite")

parser.add_option("--eps_s", type="float", default=1e-8, dest="eps_s")
parser.add_option("--toleranceH", type="float", default=0.1, dest="toleranceH")  # m/yr
parser.add_option("--toleranceXg", type="float", default=1e-3, dest="toleranceXg")  # m/yr
parser.add_option("--toleranceInner", type="float", default=1e-3, dest="toleranceInner")
parser.add_option("--initUTolerance", type="float", default=1e-3, dest="initUTolerance")
parser.add_option("--dt", type="float", default=3e-4, dest="dt")
parser.add_option("--xgInit", type="float", default=1.0, dest="xgInit")

parser.add_option("--plot", action="store_true", dest="plot")
parser.add_option("--plotContinuous", action="store_true", default=False, dest="plotContinuous")
parser.add_option("--disableLongi", action="store_true", default=False, dest="disableLongi")
parser.add_option("--useChannel", action="store_true",  default=False, dest="useChannel")
parser.add_option("--useSchoofBasal", action="store_true",  default=False, dest="useSchoofBasal")

options, args = parser.parse_args()

stepsPerWrite = options.stepsPerWrite
maxSteps = options.maxSteps

solver = Solver(options)

if(not os.path.exists(options.folder)):
  os.mkdir(options.folder)

if(maxSteps == 0):
  exit()

solver.runTimeStepping(maxSteps,stepsPerWrite)




