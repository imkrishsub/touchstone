#!/usr/bin/python
import numpy
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy.linalg as linalg

from SteadyStateSolver import SteadyStateSolver

def computeBLinear(x,solver):
  shoofx = 750000.
  slope = solver.linearSlope/solver.HBar
  pente = 778.5/solver.HBar
  intercept = -720/solver.HBar

  bx = solver.linearSlope*solver.xBar/shoofx/solver.HBar
  b0 = slope*intercept/pente

  eps_b = 1e-8
  absx = numpy.sqrt(x**2 + eps_b**2)
  b = b0 + absx*bx
 
  return b
  
  
def readResults(fileName, solver):
  print "reading from: ", fileName
  filePointer = open(fileName,'rb')
  Nx = numpy.fromfile(filePointer, dtype=int, count=1)[0]
  if(Nx != solver.Nx):
    print "Bad Nx in file."
    exit(1)
  
  solver.xg = numpy.fromfile(filePointer, dtype=float, count=1)[0]
  solver.time = numpy.fromfile(filePointer, dtype=float, count=1)[0]
  solver.H = numpy.fromfile(filePointer, dtype=float, count=Nx)
  solver.u = numpy.fromfile(filePointer, dtype=float, count=Nx)
  filePointer.close()

parser = OptionParser()

parser.add_option("--p", type="float", default=0.0, dest="p")
parser.add_option("--A", type="float", default=1e-25, dest="A")
parser.add_option("--C", type="float", default=7.624e6, dest="C")
parser.add_option("--linearSlope", type="float", default=778.5, dest="linearSlope") #drop in m per 750 km, as in Schoof 2007

parser.add_option("--inFile", type="string", default="none", dest="inFile")

parser.add_option("--Nx", type="int", default=1025, dest="Nx")
parser.add_option("--eps_s", type="float", default=1e-8, dest="eps_s")

options, args = parser.parse_args()


Nx = options.Nx



solver = SteadyStateSolver(Nx)
solver.eps_s = options.eps_s
readResults(options.inFile, solver)

xg = solver.xg
x = solver.cheb.x*xg
Dx = solver.cheb.Dx/xg

solver.p = options.p
solver.A = options.A
solver.C = options.C
solver.linearSlope = options.linearSlope

solver.computeNondimConstants()

a = solver.a
n = solver.n

epsilon = solver.epsilon
Kappa = solver.Kappa
gamma = solver.gamma

delta = solver.delta

H = solver.H
b = computeBLinear(x,solver)

Hf = solver.maxS(b)/(1-delta)

bx = numpy.dot(Dx,b)
bxx = numpy.dot(Dx,bx)

Hx = numpy.dot(Dx,H)
Hxx = numpy.dot(Dx,Hx)
Hxxx = numpy.dot(Dx,Hxx)

u = solver.u
ux = numpy.dot(Dx,u)


H0 = H[0]




HFit_0 = H


for index in range(20):
  Np = HFit_0*(1 - Hf/HFit_0)**solver.p
  uFit_0 = a*x/HFit_0
  
  xIndexMatch = numpy.nonzero(Kappa*uFit_0[2:] > Np[2:]**n)[0][0]
  xt = xg-x[xIndexMatch]
  xIndexMatch = numpy.nonzero(xg-x < 2*xt)[0][0]
  print xIndexMatch, x[xIndexMatch]
  
  #abs_u = solver.absS(uFit_0)
  #Np = HFit_0*(solver.maxS(1 - Hf/HFit_0))**solver.p
  #tauB_xg = -gamma*Np*(abs_u/(Kappa*abs_u + Np**n))**(1./n)*uFit_0/abs_u
  tauB_xg = -gamma*(uFit_0)**(1./n)
  # assume tauB = -tauD
  HxFit_0 = tauB_xg/HFit_0 + bx
  M = Dx
  M[0,:] = 0
  M[0,xIndexMatch] = 1
  rhs = HxFit_0
  rhs[0] = H[xIndexMatch]#H0
  HFit_0 = linalg.solve(M,rhs)
  
xMin = -100


xMatch = x[xIndexMatch]
xt = xg-xMatch
s = (x-xMatch)/xt

NpFit_0 = Np

Np = H*(1 - Hf/H)**solver.p


H_xg = Hf
u_xg = a*x/H_xg

ux_xg = (delta/(8*epsilon)*H_xg)**n
Hx_xg = (a - ux_xg*H_xg)/u_xg

# f = a*s**3 + b*s**2 + c*s + d
# f(0) = d = f0
# f(1) = a+b+c+d = f1
# f'(0) = c = f0'
# f'(1) = 3*a+2*b+c = f1'
# 0 0 0 1
# 1 1 1 1
# 0 0 1 0
# 3 2 1 0
A = [[0, 0, 0, 1],
     [1, 1, 1, 1],
     [0, 0, 1, 0],
     [3, 2, 1, 0]]

b = [H[xIndexMatch],H_xg[-1],Hx[xIndexMatch]*xt,Hx_xg[-1]*xt]
c = linalg.solve(A,b)

HSpline = c[0]*s**3 + c[1]*s**2 + c[2]*s + c[3]

plt.figure(1)
plt.plot(x,H,'b', x, HFit_0, 'g', x, H_xg, 'k', x, HSpline, 'm')
plt.figure(2)
plt.plot(x[xMin:],H[xMin:],'b', x[xMin:], HFit_0[xMin:], 'g', x[xMin:], H_xg[xMin:], 'k', x[xMin:], HSpline[xMin:], 'm')
plt.figure(3)
plt.plot(x[xMin:],u[xMin:],'b', x[xMin:],uFit_0[xMin:],'g', x[xMin:], Np[xMin:]**n/Kappa, 'r', x[xMin:], NpFit_0[xMin:]**n/Kappa, 'm', x[xMin:], u_xg[xMin:], 'k')
plt.figure(4)
plt.plot(x[xMin:],Hx[xMin:],'b', x[xMin:], HxFit_0[xMin:], 'g', x[xMin:], Hx_xg[xMin:], 'k')
plt.show()
exit()


H_xg = solver.maxS(b[-1])/(1-delta)
u_xg = a*xg/H_xg

ux_xg = (delta/(8*epsilon)*H_xg)**n
Hx_xg = (a - ux_xg*H_xg)/u_xg

print H[-1], H_xg
print Hx[-1], Hx_xg
print u[-1], u_xg
print ux[-1], ux_xg

xt = xg-x[xIndexMatch]
xFit = xt*solver.cheb.x + (xg-xt)
DxFit = solver.cheb.Dx/xt

HFit = H[-1] - (H[xIndexMatch]-H[-1])*(xFit-xg)/xt
uFit = a*xFit/HFit

bFit = computeBLinear(xFit,solver)

HfFit = solver.maxS(bFit)/(1-delta)
for index in range(20):
  Np = HFit*(1 - HfFit/HFit)**solver.p
  HxFit = numpy.dot(DxFit,HFit)
  uxFit = numpy.dot(DxFit,uFit)

  tauD = -HFit*(HxFit - bx[-1])
  tauB = -(gamma/Kappa**(1/n))*Np
  tauL = numpy.dot(DxFit, 4*epsilon*HFit*uxFit**(1/n))
  
  M = numpy.dot(Dx, numpy.dot(numpy.diag(4*epsilon*HFit*uxFit**(1/n-1)),Dx))
  M[0,:] = 0
  M[0,1] = 1
  M[-1,:] = Dx[-1,:]
  rhs = -tauB - tauD
  rhs[0] = 0
  rhs[-1] = ux_xg
  uFit = linalg.solve(M,rhs)
  HFit = a*x/uFit
  uxFit = numpy.dot(Dx,uFit)
  
  plt.figure(1)
  plt.plot(x,H,'b', xFit, HFit, 'g')
  plt.xlim([xFit[0],xFit[-1]])
  plt.figure(2)
  plt.plot(x,u,'b', xFit, uFit, 'g')
  plt.xlim([xFit[0],xFit[-1]])
  plt.figure(3)
  plt.plot(x,ux,'b', xFit, uxFit, 'g')
  plt.xlim([xFit[0],xFit[-1]])
  plt.figure(4)
  plt.plot(xFit,tauD,'g', xFit, tauB, 'r', xFit, tauL,'b')
  plt.show()



xMin = -100
plt.figure(1)
plt.plot(x,H,'b', x, HFit_0, 'g')
plt.figure(2)
plt.plot(x[xMin:],H[xMin:],'b', x[xMin:], HFit_0[xMin:], 'g')
plt.figure(3)
plt.plot(x[xMin:],Kappa*uFit_0[xMin:],'b', x[xMin:], Np[xMin:]**n, 'g')
plt.show()

exit()


  
H_xg = solver.maxS(b[-1])/(1-delta)
u_xg = a*xg/H_xg

ux_xg = (delta/(8*epsilon)*H_xg)**n
Hx_xg = (a - ux_xg*H_xg)/u_xg

tauD_xg = -H_xg*(Hx_xg - bx[-1])

abs_u = solver.absS(u_xg)
Np = H_xg*(solver.maxS(0.0))**solver.p
tauB_xg = -gamma*Np*(abs_u/(Kappa*abs_u + Np**n))**(1./n)*u_xg/abs_u

print tauB_xg, -gamma*Np/Kappa**(1/n)

tauL_xg = -tauB_xg - tauD_xg

uxx_xg = (tauL_xg/(4*epsilon) - Hx_xg*ux_xg**(1/n))*n*ux_xg**(1-1/n)/H_xg

Hxx_xg = (-uxx_xg*H_xg - 2*ux_xg*Hx_xg)/u_xg

p = solver.p
tauDx_xg = -Hx_xg*(Hx_xg - bx[-1]) - H_xg*(Hxx_xg - bxx[-1])
Hfx_xg = bx[-1]/(1-delta)
if((p == 0) or (p == 1)):
  if(p == 0):
    tauBx_xg = -gamma*abs_u**(1/n-1)*ux_xg
  else:
    Npx = Hx_xg*(solver.maxS(0.0))**solver.p \
      + solver.p*(solver.maxS(0.0))**(solver.p-1)*(Hx_xg-Hfx_xg)
    #tauBx_xg = -tauB_xg/n*(n*Kappa*abs_u**2*Np**(n-1)*Npx + ux_xg*Np**(2*n))/(abs_u*Np**n*(Kappa*abs_u + Np**n))
    tauBx_xg = -gamma*Npx/Kappa**(1/n)

  tauLx_xg = -tauBx_xg - tauDx_xg
  
  uxxx_xg = (tauLx_xg/(4*epsilon) - Hxx_xg*ux_xg**(1/n) - 2*Hx_xg*ux_xg**(1/n-1)*uxx_xg
             - H_xg/n*(1/n-1)*ux_xg**(1/n-2)*uxx_xg**2)*n*ux_xg**(1-1/n)/H_xg
             
  Hxxx_xg = (-uxxx_xg*H_xg - 3*uxx_xg*Hx_xg - 3*ux_xg*Hxx_xg)/u_xg

  xPrime = x-xg
  HFit_xg = H_xg + Hx_xg*xPrime + (1/2)*Hxx_xg*xPrime**2 + (1/6)*Hxxx_xg*xPrime**3
  HxFit_xg = Hx_xg + Hxx_xg*xPrime + (1/2)*Hxxx_xg*xPrime**2
  HxxFit_xg = Hxx_xg + Hxxx_xg*xPrime
  HxxxFit_xg = Hxxx_xg*numpy.ones(x.shape)

else:
  Npxp = H_xg*((Hfx_xg-Hx_xg)/H_xg)**p # xp = (xg-x)**p
  
  tauBxp_xg = -gamma*Npxp/Kappa**(1/n)
  tauLxp_xg = -tauBxp_xg

  uxxxp_xg = (-tauBxp_xg/(4*epsilon))*n*ux_xg**(1-1/n)/H_xg
  Hxxxp_xg = (-uxxxp_xg*H_xg)/u_xg

  xPrime = x-xg
  HFit_xg = H_xg + Hx_xg*xPrime + (1/2)*Hxx_xg*xPrime**2 + (1/(2+p))*(1/(1+p))*Hxxxp_xg*(-xPrime)**(2+p)
  HxFit_xg = Hx_xg + Hxx_xg*xPrime + (1/2)*Hxxx_xg*xPrime**2
  HxxFit_xg = Hxx_xg + Hxxx_xg*xPrime
  HxxxFit_xg = Hxxx_xg*numpy.ones(x.shape)



HxFit_0 = numpy.dot(Dx,HFit_0)
HxxFit_0 = numpy.dot(Dx,HxFit_0)
HxxxFit_0 = numpy.dot(Dx,HxxFit_0)

#u = solver.u
Np = HFit_0*(1 - Hf/HFit_0)**solver.p

chi = Kappa*uFit_0/Np**n
alpha1 = (1/(chi + 1))**(1/n)
alpha2 = (1/(chi + 1))
alpha3 = (Hxx-HxxFit_xg)/(HxxFit_0-HxxFit_xg)
Hxx_interp1 = alpha1*HxxFit_0 + (1-alpha1)*HxxFit_xg
Hxx_interp2 = alpha2*HxxFit_0 + (1-alpha2)*HxxFit_xg


xMin = -100

plt.figure(1)
plt.plot(x,H,'b', x, HFit_0, 'g', x, HFit_xg, 'r')
plt.figure(2)
plt.plot(x[xMin:],H[xMin:],'b', x[xMin:], HFit_0[xMin:], 'g', x[xMin:], HFit_xg[xMin:], 'r')
plt.figure(3)
plt.plot(x[xMin:],Hx[xMin:],'b', x[xMin:], HxFit_0[xMin:], 'g', x[xMin:], HxFit_xg[xMin:], 'r')
plt.figure(4)
plt.plot(x[xMin:],Hxx[xMin:],'b', x[xMin:], HxxFit_0[xMin:], 'g', x[xMin:], HxxFit_xg[xMin:], 'r',
         x[xMin:],Hxx_interp1[xMin:],'k',x[xMin:],Hxx_interp2[xMin:],'m')
         
plt.figure(5)
plt.plot(x,solver.u,'b', x, a*x/H, 'k', x, a*x/HFit_0, 'g', x, a*x/HFit_xg, 'r')
#plt.plot(x,u*Hx + ux*H,'b', x, a*numpy.ones(x.shape), 'k')

#plt.figure(5)
#plt.plot(x[xMin:],Hxxx[xMin:],'b', x[xMin:], HxxxFit_0[xMin:], 'g', x[xMin:], HxxxFit_xg[xMin:], 'r')
plt.figure(6)
plt.plot(x[xMin:],Kappa*u[xMin:],'b', x[xMin:], Np[xMin:]**n, 'g')
#plt.plot(x[xMin:],alpha1[xMin:],'k',x[xMin:],alpha2[xMin:],'m', x[xMin:],alpha3[xMin:],'b')
plt.show()

#
#H0 = H[0]
#u0 = 0
#Hx0 = -gamma*u0**(1/n)/H0 + bx[0]
#ux0 = (a - u0*Hx0)/H0
#Hxx0 = -n*gamma/H*(a/H0)**(1/n)
#HFit_0 = H0 + Hx0*x + 1/(1/n+1)*Hxx0*x**(1/n+1)
#HxFit_0 = Hx0 + Hxx0*x**(1/n)
#HxxFit_0 = Hxx0/n*x**(1/n-1)
#
#plt.figure(1)
#plt.plot(x,H,'b', x, HFit_xg, 'r', x, HFit_0, 'g')
#plt.ylim([numpy.amin(H),numpy.amax(H)])
#plt.figure(2)
#plt.plot(x,Hx,'b', x, HxFit_xg, 'r', x, HxFit_0, 'g')
#plt.ylim([numpy.amin(Hx),numpy.amax(Hx)])
#plt.figure(3)
#plt.plot(x,Hxx,'b', x, HxxFit_xg, 'r', x, HxxFit_0, 'g')
#plt.ylim([numpy.amin(Hxx),numpy.amax(Hxx)])
#plt.show()