#!/usr/bin/python
import numpy
import scipy.optimize

from scipy.integrate import odeint  
  
def bcRes(xg,solver):
  
  options = solver.options
  (b,bx) = solver.computeB(xg)
  H = b/(1-solver.delta)
  u = options.a*xg/H
  ux = (solver.delta/(8*solver.epsilon)*H)**options.n
  Hx = (options.a - ux*H)/u
  res = -solver.gamma*u**(1/options.n)-H*(Hx-bx)
  return res


def dy_dx(y, x, solver):
  options = solver.options
  H = y[0]
  u = options.a*x/H
  
  (b,bx) = solver.computeB(x)
  tauB = -solver.gamma*numpy.abs(u)**(1/options.n-1)*u
  Hx = tauB/H + bx
  return [Hx]


def initHBounded(solver, xgMin, xgMax):

  xg = scipy.optimize.brentq(bcRes, xgMin, xgMax,args=(solver,))
  return initH(solver, xg)

def initH(solver, xg):
  
  (b_xg,bx_xg) = solver.computeB(xg)
  
  # initial conditions:
  H_xg = b_xg/(1-solver.delta)
  
  y_xg = [H_xg]
  
  xInteg = solver.cheb.x[::-1]*xg
  
  y = odeint(dy_dx, y_xg, xInteg, args=(solver,))

  H = y[::-1,0]
  return (H,xg)
