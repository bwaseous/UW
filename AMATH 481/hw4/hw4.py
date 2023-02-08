import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time

from scipy import sparse
from scipy import optimize
from scipy.sparse import linalg

## 1 ##
alpha = 2
L = 10
m = 501
n = 128
xspan = np.linspace(-L,L,n+1)[:-1]
tspan = np.linspace(0,2,m)

dt = tspan[1] - tspan[0]
dx = xspan[1] - xspan[0]

u_0 = 10*np.cos(2*np.pi*xspan/L) + 30*np.cos(8*np.pi*xspan/L)
lmbd_star = alpha*dt/(dx**2)

####
def g_14(z,lmbd):
    g = 1 + (lmbd/6)*(16*np.cos(z) - np.cos(2*z) - 15)
    return g

A1 = np.abs(g_14(1,lmbd_star))

onefour_star_g = lambda z: g_14(z,lmbd_star)
z_14_opt= sp.optimize.fminbound(lambda z: -onefour_star_g(z), -np.pi, np.pi)

A2 = np.abs(g_14(z_14_opt, lmbd_star))

#####
e1 = np.ones(n)
D4 = sp.sparse.spdiags([16*e1, -1*e1,-1*e1,16*e1,-30*e1,16*e1, -1*e1,-1*e1,16*e1],
                       [-n+1,-n+2,-2,-1,0,1,2,n-2,n-1], n,n, format ='csc')/12

A3 = D4.todense()

####
u_14 = np.zeros((len(xspan),len(tspan)))
u_14[:,0] = u_0

tic = time.time()
for i in range(len(tspan)-1):
    u_14[:,i+1] = u_14[:,i] + lmbd_star*D4@u_14[:,i]
time_14 = time.time() - tic
    
A5 = u_14[:,-1].reshape(-1,1)

## 2 ##
def g_cn(z,lmbd):
    g = (lmbd*np.cos(z) + (1-lmbd))/((1+lmbd) - lmbd*np.cos(z))
    return g

cn_star_g = lambda z: g_cn(z,lmbd_star)
z_cn_opt= sp.optimize.fminbound(lambda z: -cn_star_g(z), -np.pi, np.pi)

A6 = np.abs(g_cn(z_cn_opt, lmbd_star))

####
B = sp.sparse.spdiags([-lmbd_star/2*e1, (1+lmbd_star)*e1, -lmbd_star/2*e1],
                      [-1,0,1],n,n,format='csc')
B[-1,0] = -lmbd_star/2
B[0,-1] = -lmbd_star/2

C = sp.sparse.spdiags([lmbd_star/2*e1, (1-lmbd_star)*e1, lmbd_star/2*e1],
                      [-1,0,1],n,n,format='csc')
C[-1,0] = lmbd_star/2
C[0,-1] = lmbd_star/2

A7 = B.todense()
A8 = C.todense()

####
u_cn = np.zeros((len(xspan), len(tspan)))
u_cn[:,0] = u_0

lu_solver = sp.sparse.linalg.splu(B)

tic = time.time()
for i in range(len(tspan)-1):
    u_cn[:,i+1] = lu_solver.solve(C@u_cn[:,i])
time_cn = time.time() - tic

A9 = u_cn[:,-1].reshape(-1,1)

####
u_bc = np.zeros((len(xspan), len(tspan)))
u_bc[:,0] = u_0

tic = time.time()
for i in range(len(tspan)-1):
    u_bc[:,i+1], _ = sp.sparse.linalg.bicgstab(B, C@u_bc[:,i], x0 = u_bc[:,i])
time_bc = time.time() - tic

A10 = u_bc[:,-1].reshape(-1,1)

drive = 'E'

## 3 ##
exact_128 = np.genfromtxt('exact_128.csv',delimiter = ',')
exact_256 = np.genfromtxt('exact_256.csv',delimiter = ',')

A11 = np.linalg.norm(exact_128 - A5.T)
A12 = np.linalg.norm(exact_128 - A9.T)

####
alpha = 2
L = 10
n = 256

xspan = np.linspace(-L,L,n+1)[:-1]
dx = xspan[1]-xspan[0]
dt = lmbd_star*(dx**2)/alpha
tspan = np.arange(0,2+dt,dt)

u_0 = 10*np.cos(2*np.pi*xspan/L) + 30*np.cos(8*np.pi*xspan/L)

e1 = np.ones(n)
D4_256 = sp.sparse.spdiags([16*e1, -1*e1,-1*e1,16*e1,-30*e1,16*e1, -1*e1,-1*e1,16*e1],
                       [-n+1,-n+2,-2,-1,0,1,2,n-2,n-1], n,n, format ='csc')/12

u_14_256 = np.zeros((len(xspan),len(tspan)))
u_14_256[:,0] = u_0

tic = time.time()
for i in range(len(tspan)-1):
    u_14_256[:,i+1] = u_14_256[:,i] + lmbd_star*D4_256@u_14_256[:,i]
time_14_256 = time.time() - tic

A13 = np.linalg.norm(exact_256 - u_14_256[:,-1])

####
B_256 = sp.sparse.spdiags([-lmbd_star/2*e1, (1+lmbd_star)*e1, -lmbd_star/2*e1],
                      [-1,0,1],n,n,format='csc')
B_256[-1,0] = -lmbd_star/2
B_256[0,-1] = -lmbd_star/2

C_256 = sp.sparse.spdiags([lmbd_star/2*e1, (1-lmbd_star)*e1, lmbd_star/2*e1],
                      [-1,0,1],n,n,format='csc')
C_256[-1,0] = lmbd_star/2
C_256[0,-1] = lmbd_star/2

u_cn_256 = np.zeros((len(xspan), len(tspan)))
u_cn_256[:,0] = u_0

lu_solver_256 = sp.sparse.linalg.splu(B_256)

tic = time.time()
for i in range(len(tspan)-1):
    u_cn_256[:,i+1] = lu_solver_256.solve(C_256@u_cn_256[:,i])
time_cn_256 = time.time() - tic

A14 = np.linalg.norm(exact_256 - u_cn_256[:,-1])