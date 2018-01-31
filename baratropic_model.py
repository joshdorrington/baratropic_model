import numpy as np
from numpy.random import RandomState
import math
import time as tm
import scipy.sparse as sps
import pandas as pd

t1=tm.time()

#these functions return the values of certain coefficients in the model
def a(m,b):
    return 8.0*np.sqrt(2)*(m^2)*(b**2 + m**2 -1)/(math.pi*(-1+ 4.0*m**2)*(b**2 + m**2))


def beta(m,b,B):
    return B*b**2/(b**2 + m**2);


def gammaprime(m,b,g):
    return g*4*np.sqrt(2)*b*m/((4*(m**2)-1)*math.pi);

              
def gamma(m,b,g):
    return g*4*np.sqrt(2)*b*(m**3)/((4*(m**2)-1)*math.pi*(b**2 + m**2));


def d(m,b):
    return 64*np.sqrt(2)*(b**2 - m**2 +1)/(15*math.pi*(b**2 + m**2));
    
    

x_init=[0.927796,0.160574,-0.0969953,-0.672758,-0.118184,0.214505] #the initial conditions of the model (these are arbitrary)
dt=2e-4 #the time step used by the euler maruyama scheme
Nsample=int(5e3) #the temporal sampling rate of the data in units of dt
Nmax=1000 #how many samples you want in the final data
sigma=0.0 #the stochastic forcing coefficient
save=True # True if you wish to save your data to file, False otherwise.

xdot = np.zeros(6);


#constant def list:
C=0.1 #Newtonian damping towards the forced streamfunction xf
B=1.25 # Beta, the coriolis terms
g=0.2 #gamma, a topographic height parameter
b=0.5 #ratio of channel length to width
e=16.0*np.sqrt(2.0)/(5.0*math.pi) #a nonlinear advection coefficient

a1=a(1,b)
a2=a(2,b)
beta1=beta(1,b,B)
beta2=beta(2,b,B)
gamma1=gamma(1,b,g)
gamma2=gamma(2,b,g)
gammaprime1=gammaprime(1,b,g)
gammaprime2=gammaprime(2,b,g)
d1=d(1,b)
d2=d(2,b)

#defines the forcing vector
xf=np.zeros_like(xdot) 
xf[0]=0.95
xf[3]=-0.76095

#initialises data structures for the various parts of the system
state=np.zeros([Nmax+1,len(xdot)])
time=np.zeros(Nmax+1)
state[0]=x_init;
x=np.zeros_like(xdot)
dW=np.zeros_like(xdot)

#builds the linear operator
Lin_Op=sps.dok_matrix((len(xdot),len(xdot)))
Lin_Op[0,:]=[-C,0.0,gammaprime1,0.0,0.0,0.0]
Lin_Op[1,:]=[0.0,-C,beta1,0.0,0.0,0.0]
Lin_Op[2,:]=[-gamma1,-beta1,-C,0.0,0.0,0.0]
Lin_Op[3,:]=[0.0,0.0,0.0,-C,0.0,gammaprime2]
Lin_Op[4,:]=[0.0,0.0,0.0,0.0,-C,beta2]
Lin_Op[5,:]=[0.0,0.0,0.0,-gamma2,-beta2,-C]
Lin_Op=Lin_Op.tocsr()

#seed random numbers in thread-safe way
myrng=RandomState(int(tm.time()))

t2=tm.time()
lineartime=0
stochastictime=0
nonlineartime=0   
for N in range(1,Nmax):
    
    x=state[N-1]
    t4=tm.time() #TIME CHECK
    dW=myrng.normal(0.0,np.sqrt(dt),[Nsample,len(xdot)])#initialise random numbers for next set of timesteps
    stochastictime+=tm.time()-t4 #TIME CHECK
    for t in range(0,Nsample):
        
        #linear part done with matrix method
        t4=tm.time()
        xdot=Lin_Op.dot(x)+C*xf
        lineartime+=tm.time()-t4 #TIME CHECK
        #nonlinear part done elementwise
        t4=tm.time()
        xdot[1]+= -a1*x[0]*x[2]-d1*x[3]*x[5]
        xdot[2]+=  a1*x[0]*x[1]+d1*x[3]*x[4]
        xdot[3]+=  e*(x[1]*x[5]-x[2]*x[4])
        xdot[4]+= -a2*x[0]*x[5]-d2*x[2]*x[3]
        xdot[5]+=  a2*x[0]*x[4]+d2*x[1]*x[3]
        nonlineartime+=tm.time()-t4 #TIME CHECK
        #tendency and stochastic forcing are added to x
        x+= xdot*dt +sigma*dW[t]
    
    
    #the value of x and t is recorded after each sample loop
    time[N]=time[N-1]+(t+1)*dt
    state[N]=x
    print "%.1f %% done" % (100.0*N/Nmax)

t3=tm.time()
#writes data to the names file
if save:
    tx=[time,state]
    filepath='\\Users\\Josh\\Documents\\dphil\\baratropic_python_implementation_test.csv'
    df=pd.DataFrame(tx)
    df.to_csv(filepath)


print "Initialisation took %f seconds, the model took %d seconds" % (t2-t1,t3-t2)
print "Linear calculations took %f seconds, Nonlinear took %f seconds, stochastic took %f seconds" % (lineartime,nonlineartime,stochastictime)

