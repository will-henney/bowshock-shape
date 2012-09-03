import scipy.integrate as sci
import numpy as np
import argparse
import matplotlib.pyplot as plt

#List of functions I'll need.


def I(mu):
    """
    Function created to solve numerically the integral \int_{mu}^{1} \sqrt{x-x^3} dx
    Arguments: mu, Type: float. Description: Is the lower limit of integral
    Output: r, Type: float. Descripion: The integral's solution
    """
    f= lambda x: np.sqrt(x-x**3)
    r,err = sci.quad(f,mu,1.)
    return r

def II(mu):
    """
    This function creates a numpy array with the solutions for a given function
    Arguments: mu, Type: numpy array. Description: The array with the cosines of angle theta
               f: Function of mu
    output: ivec: type: numpy array. Description: The solutions of the given function for each element of the mu array.
    """
    ivec = np.array([I(m) for m in mu])
    return ivec

def F(mu, beta):
    """
    Function of mu, with a beta parameter. Forms part of equation:t1*cot(t1) -1 = F(t) and mu = cos(t)
    """
    R = 6.*beta*II(mu)*(mu/np.sqrt(1.-mu**2))-2.4*beta*(1.-mu**2.5)
    return R

def G(mu,beta):
    """
    The equivalent of F(mu) for the spherical case
    """
    t = np.arccos(mu)
    T = beta*(t*np.cos(t)/np.sin(t)-1.)
    return T

def Mu1(mu,beta,f):
    """
    This function solves for mu1 from the equation t1*cot(t1) -1 = F(t) 
    and approximating t1*cot(t1) ~ 1 - (t1^2)/3 - (t1^4)/45 
    with mu1 = cos(t1) 
    """
    t1=np.sqrt( 7.5*( -1.+np.sqrt( 1.-0.8*f ) ) )
    mu1=np.cos(np.pi*t1/180.)
    return mu1

#Adding external parameters
parser = argparse.ArgumentParser(description="Input parameters")
parser.add_argument("--R0D",'-R',default = 0.1, type= float, help = "R0/D ratio")
parser.add_argument("--dist",'-D',default = 5, type= float, help = "distance between proplyd and masive star")
cmdargs = parser.parse_args()

R0D = cmdargs.R0D
beta = (R0D/(1.-R0D))**2
D = cmdargs.dist

#Creating the array of mu and calculating mu1 for the spherical case and the non spherical case
#For negative mu, mu1 has a constant value. Mu has an asymtotic value. That value is such that theta + theta1 = pi with mu=cos(theta)
#mu1c= Mu1(0,beta,F(0,beta))
#print mu1c 
Mu = np.arange(0,1.05,0.05)
#Mu1 = Mu1(Mu,beta,F(Mu,beta))
Mu1Sph = Mu1(Mu,beta,G(Mu,beta))
print Mu1Sph
nanFinder = np.isnan(Mu1Sph)
#Mu1[nanFinder] = 1.  # Replace nan values for the correct value. Please ignore warnings while running the program.
Mu1Sph[nanFinder] = 1. # In both cases (spherical and non spherical), the value of mu1 is 1 when mu = 1. Check the Figure 1 from Canto,Raga & Wilkin 1996



#Calculating the radius of bowshocks in function of mu and mu1 and replacing the nan values for the correct ones (in mu =1, R is R0)

#RD = np.sqrt(1.-Mu1**2)/(Mu*np.sqrt(1.-Mu1**2)+Mu1*np.sqrt(1.-Mu**2))
RDSph = np.sqrt(1-Mu1Sph**2)/(Mu*np.sqrt(1-Mu1Sph**2)+Mu1Sph*np.sqrt(1-Mu**2))
#nanFinder = np.isnan(RD)
#RD[nanFinder] = R0D
nanFinder = np.isnan(RDSph)
RDSph[nanFinder] = R0D

print RDSph



