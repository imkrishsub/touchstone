#!/usr/bin/python
import numpy
import scipy.optimize

from scipy.integrate import odeint  
  
def bcRes(xg,solver):
  
  (b,bx) = solver.computeB(xg,solver)
  H = b/(1-solver.delta)
  u = solver.a*xg/H
  ux = (solver.delta/(8*solver.epsilon)*H)**solver.n
  Hx = (solver.a - ux*H)/u
  res = -solver.gamma*u**(1/solver.n)-H*(Hx-bx)
  return res


def dy_dx(y, x, solver):
  H = y[0]
  u = solver.a*x/H
  
  (b,bx) = solver.computeB(x,solver)
  tauB = -solver.gamma*numpy.abs(u)**(1/solver.n-1)*u
  Hx = tauB/H + bx
  return [Hx]


def initH(solver, xgMin, xgMax):

  xg = scipy.optimize.brentq(bcRes, xgMin, xgMax,args=(solver,))
  
  (b_xg,bx_xg) = solver.computeB(xg,solver)
  
  # initial conditions:
  H_xg = b_xg/(1-solver.delta)
  
  y_xg = [H_xg]
  
  xInteg = solver.cheb.x[::-1]*xg
  
  y = odeint(dy_dx, y_xg, xInteg, args=(solver,))

  H = y[::-1,0]
  return (H,xg)
