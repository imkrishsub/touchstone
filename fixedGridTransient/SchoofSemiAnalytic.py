import numpy

def computeBPoly(x):
  schoofx = 750000.
  HBar = 1000. # m hight scaling factor
  xBar = 1000000. # m domain scaling factor
  xs = x*xBar/schoofx
  b = -(729 - 2184.8*xs**2 + 1031.72*xs**4 - 151.72*xs**6)/HBar  
  return b

def computeBxPoly(x):
  schoofx = 750000.
  HBar = 1000. # m hight scaling factor
  xBar = 1000000. # m domain scaling factor
  xs = x*xBar/schoofx
  bx = -xBar/schoofx*(- 2.0*2184.8*xs + 4.0*1031.72*xs**3 - 6.0*151.72*xs**5)/HBar  
  return bx

def schoofSemiAnalyticA(xg):
  rho_i = 900.0 # kg/m^3 ice density
  rho_w = 1000.0 # kg/m^3 water density
  delta = 1.0 - rho_i/rho_w
  C = numpy.array(7.624e6)
  a = numpy.array(1.0)
  
  sPerY = 365.25*24.*3600. # number of seconds per year
  g = 9.8 # m/s^2 gravity acceleration
  aBar = .3/sPerY # m.s-1 accumulation rate
  HBar = 1000. # m hight scaling factor
  xBar = 1000000. # m domain scaling factor
  n = 3. # Glen's flow parameter
  m = 1./n
    
  q = aBar*xBar * a*xg
  
  b = HBar * computeBPoly(xg)
  bx = HBar/xBar * computeBxPoly(xg)
  H = b/(1-delta)
  u = q/H

  # -a + q/h*(bx - C/(rho_i*g)*q**m/h**(m+1)) + A/4**n * (rho_i*g*delta)**n*h**(n+1) == 0
  num = -aBar*a + u*(bx - C/(rho_i*g*H)*u**m)
  denom = (rho_i*g*delta/4.)**n*H**(n+1)
  A = -num/denom


#  plt.figure(1)
#  plt.semilogx(1/A,xBar*xgs)
#  plt.show()
#  crash
#  
  return A
  
def schoofSemiAnalyticA2(xg):
  rho_i = 900.0 # kg/m^3 ice density
  rho_w = 1000.0 # kg/m^3 water density
  delta = 1.0 - rho_i/rho_w
  C = 7.624e6
  
  sPerY = 365.25*24.*3600. # number of seconds per year
  g = 9.8 # m/s^2 gravity acceleration
  aBar = .3/sPerY # m.s-1 accumulation rate
  HBar = 1000. # m hight scaling factor
  xBar = 1000000. # m domain scaling factor
  n = 3. # Glen's flow parameter
  m = 1./n
  a = aBar
  
  q = a*xBar*xg
  
  b = HBar * computeBPoly(xg)
  bx = HBar/xBar * computeBxPoly(xg)
  H = b/(1-delta)
  u = q/H
  taub = -C*u**m
  taud = -taub 
  sx = -taud/(rho_i*g*H) # taud = -rho_i*g*H*sx
  Hx = sx+bx # sx = Hx-bx
  
  ux = (a - u*Hx)/H # H*ux + u*Hx = a
  # 2*A**(-1./n)*ux**(1./n) = 0.5*rho_i*g*delta*H
  A = ux*(0.25*rho_i*g*delta*H)**(-n)
  
  # a + q/h*(bx - C/(rho_i*g)*q**m/h**(m+1)) + A/4**n * (rho_i*g*delta)**n*h**(n+1) == 0
  #num = aBar*a + q/h*(bx - C/(rho_i*g)*q**m/h**(m+1))
  #denom = (rho_i*g*delta/4.)**n*h**(n+1)
  #A = -num/denom
#  plt.figure(1)
#  plt.semilogx(1/A,xBar*xgs)
#  plt.show()
#  crash
#  
  return A
  
def schoofSemiAnalyticB(xg):
  rho_i = 900.0 # kg/m^3 ice density
  rho_w = 1000.0 # kg/m^3 water density
  delta = 1.0 - rho_i/rho_w
  C = numpy.array(7.624e6)
  a = numpy.array(1.0)
  
  sPerY = 365.25*24.*3600. # number of seconds per year
  g = 9.8 # m/s^2 gravity acceleration
  aBar = .3/sPerY # m.s-1 accumulation rate
  HBar = 1000. # m hight scaling factor
  xBar = 1000000. # m domain scaling factor
  n = 3. # Glen's flow parameter
  m = 1./n
    
  q = aBar*xBar * a*xg
  
  b = HBar * computeBPoly(xg)
  #bx = HBar/xBar * computeBxPoly(xg)
  h = b/(1-delta)
  # q**(m+1)/h**(m+2)*(-C/(rho_i*g)) + A/4**n * (rho_i*g*delta)**n*h**(n+1) == 0
  num = -q**(m+1)/h**(m+2)*(C/(rho_i*g))
  denom = (rho_i*g*delta/4.)**n*h**(n+1)
  A = -num/denom
#  plt.figure(1)
#  plt.semilogx(1/A,xBar*xgs)
#  plt.show()
#  crash
#  
  return A

def schoofSemiAnalyticB2(xg):
  rho_i = 900.0 # kg/m^3 ice density
  rho_w = 1000.0 # kg/m^3 water density
  delta = 1.0 - rho_i/rho_w
  C = numpy.array(7.624e6)
  a = numpy.array(1.0)
  
  sPerY = 365.25*24.*3600. # number of seconds per year
  g = 9.8 # m/s^2 gravity acceleration
  aBar = .3/sPerY # m.s-1 accumulation rate
  HBar = 1000. # m hight scaling factor
  xBar = 1000000. # m domain scaling factor
  n = 3. # Glen's flow parameter
  m = 1./n
    
  q = aBar*xBar * a*xg
  
  
  b = HBar * computeBPoly(xg)
  #bx = HBar/xBar * computeBxPoly(xg)
  h = b/(1-delta)
  # q = (A*(rho_i*g)**(n+1)*delta**n/4**n/C)**(1/(m+1))*h**((m+n+3)/(m+1))
  
  A = q**(m+1)/h**(m+n+3)*4**n*C/(rho_i*g)**(n+1)/delta**n

  return A
  