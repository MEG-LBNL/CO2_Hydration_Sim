#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 2019
@author: pcagbo
"""
from scipy.integrate import odeint
import numpy as np
from tqdm import tqdm

#Empirical carbonate speciation rate/equilibrium constants.
kd = 1e10; kdp = kd/2e-3 #CO2(g) <==> CO2(aq)
k4 = 1e7; k4p = 5e10 #H2CO3 <==> proton + HCO3-
k5 = 3; k5p = 5e10 #HCO3- <==> proton + CO32-
k6 = 4e-4; k6p = 1.21e4 #HCO3- <==> CO2 + OH-
k7 = 2 #k7=diffusion-limited M(CO3) complexation; k[M2+] ~ 10 s-1
k8 = 2 #k8=diffusion-limited M(OH)2 formation.
khydrate = 0.04; kdehydrate = 12 #CO2 + H20 <==> H2CO3, uncatalyzed reaction

#Carbonic Anhydrase parameters
kcat = 4.4e6; kcatp = kcat/10 #originally kcat/2.5e-3; altering this significantly changes enzyme pH dependence
km = 1.25e-2
k1 = 1e9; k1p = km*k1#approximates that kf >> kcat.
Et = 1e-8 #CA concentration, molar.

#pH-independent solubility constants for metal-carbonates and metal-hydroxides:
Ksp1 = 3.3e-6 #Solubility product for calcium carbonate (calcite); 6.0e-9 for aragonite.
Ksp2 = 5.5e-6 #Solubility product for calcium hydroxide.
k7p = k7*Ksp1
k8p = k8*Ksp2

#Assumes a CA mechanism: CO2 + E (k1, k1p) <===> [E-CO2]* + H20 <===> E + H2CO3 (kcat, kcatp)
def HCAII(CO2aq, H2CO3):
    f = k1*CO2aq/(kcat + k1p) * (kcat/k1p + 1)
    k2 = kcat*(1 - 1/(1 + f)); k2p = kcatp * H2CO3 / (1 + f)
    return [k2, k2p, k2-k2p]

#Standard Michaelis-Menten rate equation.
def MM(CO2aq):
    return (kcat*CO2aq)/(km + CO2aq)

def ODE_system(state, t, kcat, km):
    #Defining system of ODEs:
    M, CO2g, CO2aq, MCO3, H2CO3, HCO3, CO3, proton, MOH = state
    
    #Definition for OH-:   
    pH = -np.log10(proton); pOH = 14-pH; OH = 10**-pOH

    #Definitions for complex variables f, k2 and k2p; derived using a reversible MM case
    k2, k2p = HCAII(CO2aq, H2CO3)[0:2]
    
    #system of differential eqs for carbon speciation and enzyme catalysis:
    dCO2g = 0; dCO2aq = 0; dproton = 0; dM = 0 #steady-state variables
    dH2CO3 = k2*Et - k2p*Et - k4*H2CO3 + k4p*proton*HCO3 + khydrate*CO2aq - kdehydrate*H2CO3
    dHCO3 = k4*H2CO3 - k4p*proton*HCO3 - k5*HCO3 + k5p*proton*CO3 - k6*HCO3 + k6p*OH*CO2aq
    dCO3 = -k5p*proton*CO3 + k5*HCO3 - k7*M*CO3 + k7p*MCO3
    dMCO3 = k7*M*CO3 - k7p*MCO3
    dMOH = k8*M*OH - k8p*MOH
    
    #print [t, pH, CO2aq, HCO3, MOH, MCO3]
    return [dM, dCO2g, dCO2aq, dMCO3, dH2CO3, dHCO3, dCO3, dproton, dMOH]

#DIC speciation for pH 9.5 used as inital state (pH 9.5 = 3.16e-5 M OH- ion.)
pHi = 9.0 ; Hi = 10**-pHi #initial pH ; [proton]
#Initial state concentrations of DIC species:
Ka1 = 4.45e-7; Ka2 = 4.69e-11; DICi = 1.8e-3 #apply equilibrium [DIC] for 470 ppm CO2 atmosphere ([DIC]i = [DIC]eq)
H2CO3i = DICi*Hi**2/(Hi**2 + Hi*Ka1 + Ka1*Ka2)
HCO3i = DICi*Hi*Ka1/(Hi**2 + Hi*Ka1 + Ka1*Ka2)
CO3i = DICi*Ka1*Ka2/(Hi**2 + Hi*Ka1 + Ka1*Ka2)

"""Initialize the concentrations of intermediates for t = 0.
init_state = [M2+, CO2g, CO2aq, MCO3, H2CO3, HCO3, CO3, proton, M(OH)2]"""
init_state = [0.01, 4.7e-4, 1e-5, 0, H2CO3i, HCO3i, CO3i, Hi, 0] #assume at T=0 system starts in non-equilibrium state.
dt = 10; t = np.arange(0, 3610, dt) #run 24 hr simulation

state = odeint(ODE_system, init_state, t, args=(kcat, km)) #state matrix of concentrations

def state_matrix():
    t2 = np.arange(0, t[-1], dt)
    return t2, odeint(ODE_system, init_state, t2, args=(kcat, km))