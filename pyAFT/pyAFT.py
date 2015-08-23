import numpy as np
from .parameters import MAXCOUNT
from .utilities import drawbinom
import pylab as plt
import random


def generate_synthetic_counts(rho, Nc=30):
    """Generate Synthetic AFT data.

    Parameters:
    rho : track density
    Nc : Number of crystals

    """
    # Probability in binomial distribution
    prob = rho / (1. + rho)

    # For Nc crystals, generate synthetic Ns and Ni count data using binomial
    # distribution, conditional on total counts Ns + Ni, sampled randomly with
    # a maximum of 1000.
    # Nc is the number of
    NsNi = np.random.randint(0, MAXCOUNT, Nc)
    Ns = np.array([drawbinom(I, prob) for I in NsNi])
    Ni = NsNi - Ns
    return {"Spontaneous tracks (Ns)": Ns, "Induced tracks (Ni)": Ni}


def calculate_central_age(Ns, Ni, zeta, seZeta, rhod, seRhod, sigma=0.15):
    """Function to calculate central age."""

    #Calculate mj
    lbda = 1.55125e-10
    m  = Ns+Ni
    p  = Ns/m
    
    theta = sum(Ns)/sum(m)
    
    
    for i in range(0, 30):
        w = m/(theta*(1-theta)+(m-1)*theta**2*(1-theta)**2*sigma**2)
        #Calculate new value of sigma and theta
        sigma = sigma*sqrt(sum(w**2*(p-theta)**2)/sum(w))
        theta = sum(w*p)/sum(w)
    
    t = (1/lbda)*log(1+1/2*lbda*zeta*rhod*(theta)/(1-theta))
    se = sqrt(1/(theta**2*(1-theta)**2*sum(w))+(seRhod/rhod)**2+(seZeta/zeta)**2)*t
    
    return {"Central":t/1e6, "2 sigma Error": 2*se/1e6, "Dispersion":sigma*100}

def calculate_pooled_age(Ns, Ni, zeta, seZeta, rhod, seRhod):

    #Calculate mj
    lbda = 1.55125e-10
    sigma = 0
    m = Ns+Ni
    p = Ns/m
    
    theta = sum(Ns)/sum(m)
    
    
    for i in range(0, 30):
        w = m/(theta*(1-theta)+(m-1)*theta**2*(1-theta)**2*sigma**2)
        theta = sum(w*p)/sum(w)
    
    t = (1/lbda)*log(1+1/2*lbda*zeta*rhod*(theta)/(1-theta))
    se<-sqrt(1/sum(Ns)+1/sum(Ni)+(seRhod/rhod)**2+(seZeta/zeta)**2)*t
    
    return {"Pooled Age":t/1e6, "2 sigma Error": 2*se/1e6}


def calculate_single_grain_ages(Ns, Ni, rhod, zeta, g=0.5, trf="Linear"):
    # Total Decay constant for 238U
    lbda = 1.55125e-10

    # Linear Transformation
    if (trf == "linear"):
        z = 1/lbda*log(1+g*zeta*lbda*rhod*(Ns/Ni))
        sez = z*np.sqrt(1/Ns + 1/Ni)
        z0 = (sum(z/sez**2) / sum(1/sez**2)

    #Logarithmic transformation
    if trf is "Log":
        z = log(g*zeta*lbda*rhod*(Ns/Ni))
        sez = z*sqrt(1/Ns+1/Ni)
        z0 = log(g*zeta*lbda*rhod*(sum(Ns)/sum(Ni)))
    
    #Arcsine
    if trf is "arcsine":
        z = asin(sqrt((Ns+3/8)/(Ns+Ni+3/4)))
        sez = 1/(2*sqrt(Ns+Ni))
        z0 = asin(sqrt(sum(Ns)/sum(Ns+Ni)))


    Age = z/1e6
    Error = sez / 1e6
    
    return Age, Error

def radialplot(Ns, Ni, zeta, seZeta, rhoD, seRhoD):

    #Calculate single grain ages
    Ages, Errors = calculate_single_grain_ages(Ns, Ni, rhoD, zeta)
    
    #Calculate central age
    CentralAge = calculate_central_age(Ns, Ni, zeta, seZeta, rhoD, seRhoD)
    #Calculate pooled age
    PooledAge = calculate_pooled_age(Ns, Ni, zeta, seZeta, rhoD, seRhoD)
    
    x = Ages / Errors
    y = (Ages - CentralAge) / Errors
    
    A = plt.plot(x, y, "o")
    return A


def PlotPlot(age, error):
     
    # scale
    z = age
    ze = error
    z0 = sum(age / error**2) / sum(1 / error**2)      



def project_fission_track()
    """Project fission track onto C-axis"""
    return



