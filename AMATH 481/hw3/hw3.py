import numpy as np
import scipy as sp
from scipy import sparse
from scipy import integrate
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import time

## Problem 1 ##
L = 10
h = 0.1
xspan = np.linspace(-L, L, 201)[:-1]
tspan = np.arange(0,10.5,.5)
u_0 = np.exp(-1*(xspan - 5)**2)

## a ##
A = sp.sparse.spdiags([-1*np.ones(len(xspan)), np.ones(len(xspan))],
                      [-1,1], len(xspan), len(xspan), format = 'csc')/(2*h)
A[-1,0] = 1/(2*h)
A[0,-1] = -1/(2*h)
A1 = A.todense()

## b ##
def advection(t,x,A,c):
    x = -c*(A@x)
    return x

c = -0.5
sol = sp.integrate.solve_ivp(advection, [0,10], args = (A, c), y0 = u_0, 
                             t_eval = tspan, method = 'RK45')
A2 = sol.y

## c ##
def heaviside_advection(t,x,A,xspan):
    c = -(1 + 2*np.sin(5*t) - np.heaviside(xspan-4,0))
    x = -c*(A@x)
    return x

sol = sp.integrate.solve_ivp(heaviside_advection, [0,10], args = (A, xspan), 
                             y0 = u_0, t_eval = tspan, method = 'RK45')
A3 = sol.y

## Problem 2 ##
nu = 0.001
L = 10
m = 64
n = m**2
dt = 0.5
t_max = 4

xspan = np.linspace(-L,L,m+1)[:-1]
yspan = np.linspace(-L,L,m+1)[:-1]
tspan = np.arange(0,t_max + dt,dt)

h = xspan[1] - xspan[0]

X, Y = np.meshgrid(xspan, yspan)
w_0 = np.exp(-2*X.T**2 - (Y.T**2)/20).flatten()

## a ##
e1 = np.ones(n) # vector of ones
Low1 = np.tile(np.concatenate((np.ones(m-1), [0])), (m,)) # Lower diagonal 1
Low2 = np.tile(np.concatenate(([1], np.zeros(m-1))), (m,)) #Lower diagonal 2

Up1 = np.roll(Low1, 1) # Shift the array for spdiags
Up2 = np.roll(Low2, m-1) # Shift the other array

A = sp.sparse.spdiags([e1, e1, Low2, Low1, -4*e1, Up1, Up2, e1, e1],
            [-(n-m), -m, -m+1, -1, 0, 1, m-1, m, (n-m)], n, n, format = 'csc')/h**2
A[0,0] = 2/h**2

B = sp.sparse.spdiags([e1, -e1, e1, -e1], [-(n - m), -m, m, (n-m)], 
                      n, n, format = 'csc')/(2*h)

C = sp.sparse.spdiags([Low2, -Low1, Up1, -Up2], [-m+1,-1,1,m-1], n, n, format = 'csc')/(2*h)

A4 = A.todense()
A5 = B.todense()
A6 = C.todense()

## bi ##
def gauss_vort(t, x, A, B, C, nu):
    psi = sp.sparse.linalg.spsolve(A,x)
    x = -((B@psi)*(C@x) - (C@psi)*(B@x)) + nu*(A@x)
    return x

t_gauss = time.time()
sol = sp.integrate.solve_ivp(gauss_vort, [0,t_max], y0 = w_0, args = (A,B,C,nu), 
                            t_eval = tspan, method = 'RK45')
t_gauss = time.time() - t_gauss

A7 = sol.y.T

## bii ##
lusolver = sp.sparse.linalg.splu(A)

def plu_vort(t, x, lusolver, B, C, nu):
    psi = lusolver.solve(x)
    x = -((B@psi)*(C@x) - (C@psi)*(B@x)) + nu*(A@x)
    return x

t_plu= time.time()
sol = sp.integrate.solve_ivp(plu_vort, [0,t_max], y0 = w_0, args = (lusolver,B,C,nu), 
                            t_eval = tspan, method = 'RK45')
t_plu = time.time() - t_plu

A8 = sol.y.T

A9 = sol.y.T.reshape(len(tspan),m,m)