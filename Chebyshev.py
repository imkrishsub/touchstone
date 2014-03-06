#!/usr/bin/python
import numpy
import math
#import scipy.fftpack
#import numpy.linalg as linalg
import chebtran


def transformPhysToCheb(field):
  return chebtran.ifcglt(field)

def transformChebToPhys(field_f):
  return chebtran.fcglt(field_f)
  
def intXCheb(field_f,xMax):
  Nx = field_f.size
  field_int = numpy.zeros(field_f.shape,float)
  field_int[-1] = field_f[-2]/(2.0*(Nx-1))
  field_int[1] = field_f[0] - 0.5*field_f[2]
  for n in range(2,Nx-1):
    field_int[n] = (-field_f[n+1]+field_f[n-1])/(2.0*n)
  field_int *= -xMax/2.0
  
  return field_int


def dxCheb(field_f,xMax):
  Nx = field_f.size
  temp = 2.0*numpy.arange(Nx)*field_f
  field_dx = numpy.zeros(field_f.shape,float)
  field_dx[Nx-3:Nx-1]=temp[Nx-2:Nx]   
  for j in numpy.arange(numpy.floor(Nx/2)-1):
    k=Nx-3-2*j
    field_dx[k] = temp[k+1] + field_dx[k+2]
    field_dx[k-1] = temp[k] + field_dx[k+1]
  field_dx[0] = (temp[1] + field_dx[2])*0.5

  field_dx *= -2.0/xMax
  return field_dx
 

def intXPhys(field,xMax):
  result = transformChebToPhys(intXCheb(transformPhysToCheb(field),xMax))
  result = result - result[0]
  return result
  
 
def dxPhys(field,xMax):
  return transformChebToPhys(dxCheb(transformPhysToCheb(field),xMax))


def antialias(field):
  Nx = field.size
  nMax = int(2*Nx/3)
  field_f = transformPhysToCheb(field)
  field_f[nMax:] = 0.
  return transformChebToPhys(field_f)
  
def null(A, eps=1e-15):
  import scipy
  from scipy import matrix
  A = matrix(A)
  u, s, vh = scipy.linalg.svd(A)
  null_mask = (s <= eps)
  null_space = scipy.compress(null_mask, vh, axis=0)
  return scipy.transpose(null_space)
  
class Chebyshev:
  def __init__(self, Nx):
    #self.Tx = transformChebToPhys(numpy.identity(Nx)) #numpy.zeros((Nx,Nx),float)
    #self.invTx = transformPhysToCheb(numpy.identity(Nx)) #numpy.zeros((Nx,Nx),float)
    self.Dx = numpy.zeros((Nx,Nx),float)
    #self.Dxx = numpy.zeros((Nx,Nx),float)
    #self.intX = numpy.zeros((Nx,Nx),float)
    for n in range(Nx):
      delta = numpy.zeros((Nx),float)
      delta[n] = 1.0
      delta_f = transformPhysToCheb(delta)
      #self.Tx[:,n] = transformChebToPhys(delta)
      #self.invTx[:,n] = delta_f
      deltax_f = dxCheb(delta_f,1.0)
      #self.intX[:,n] = intXPhys(delta,1.0)
      self.Dx[:,n] = transformChebToPhys(deltax_f)
      #self.Dxx[:,n] = transformChebToPhys(dxCheb(deltax_f,1.0))
    
#    print self.intX[-10:,-10:]
#    
#    v = null(self.intX)
#    print v
#    v = null(self.intX[0:-1,1:])
#    print v
#    inv = linalg.inv(self.intX[0:-1,1:])
#    self.intX[:,:] = 0.0
#    self.intX[1:,0:-1] = inv
#    
#    print self.intX[0:10,0:10]
#    self.intX = numpy.dot(self.invTx,numpy.dot(self.intX ,self.Tx))
    self.x = 0.5*(1 - numpy.cos(math.pi*numpy.arange(Nx,dtype=float)/float(Nx-1)))

    