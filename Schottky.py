# Code to Simulate a lateral p-type Schottky diode device
#----------------------------------------------------------
#Imports
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
#----------------------------------------------------------
#CONSTANTS
#q is the charge on an electron (C)
q=1.60E-19
#me is the mass of an electron (kg)
me=9.11E-31
#h is planck's constant (kg m^2 s^-2 s)
h= 6.63E-34
#k is the Boltzman costant (m^2 kg s^-2 K^-1 or JK^-1)
k=1.38E-23

#Material Parameters
#METAL
#wf_m is the metal work function (eV)
wf_m=3.8

#SEMICONDUCTOR
#Semiconductor Bandgap(eV)
Eg=5.47
#wf_s is the work function of the semiconductor(eV)
wf_s=4.2
#X_s is the electron affinity of the semiconductor(eV)
X_s=1
#Na is the doping concentration m^-3
Na=1.2E24
#p_s is the permitivity of the semiconductor (m^-3 kg^-1 s^4 A^2)
p_s=5.51
#p_0 is the permitivity of free space (m^-3kg^-1s^4A^2)
p_0=8.85E-12
#hole mobility in boron doped diamond (m^2 V^-1 s^-1)
mu_p=0.2
#Na is the p-type semiconductor acceptor density
N_a=1.2E24
#Internal reistance of the semiconductor (ohms)
R_int=4.85E-10
#A is the diode area m^2 (the size of the schottky contact pi*90nm^2 in m^2)
A=2.54E-14
#L is the schottky diode channel length in m
L= 20E-9
#R is the effective Richardson constant A m^-2 K^-2 given in [1]
R=90.0E4
#T is the absolute temperature i.e. room temperature (K)
T=296.0 
#n is the ideality factor
n=1.008

#----------------------------------------------------
#CALCULATIONS

#SB_0 is the Schottky barrier at zero bias, this value has been reported in [4]
#SB_0=0.196 
#SB_0 is the Schottky barrier at zero bias [2]
SB_0=Eg-(wf_m-X_s)

#The following equations determine the saturation current: [1],[2],[3]   
#exponential for saturation current
expon=math.exp((-q*SB_0)/(k*T))
#I_0 is the saturation current which is related to the schottky barrier height at zero bias
I_0=R*T**2*expon     


#----------------------------------------------------------
#DETERMINATION OF THE CURRENT


#V_threshold is the voltage at which the diffusion current dominates
V_t=wf_s-wf_m
print("The threshold voltage is")
print(V_t)

#V_er is the voltage after which ohmic transport dominates [5]

V_er=(4/3)*((q*N_a*L**2)/(p_s*p_0))
print("V_er is")
print(V_er)

#m is in the Mott-Gurney current equation and is usually between 2.3 and 3
m=3


V=0.0
I=0.0
V_arr=[]
I_arr=[]


for i in range (0,1001):
    
    V=i/100.0

    if V<=V_t:
 #I is the current found in [1],[2],[3]  for the first term and [5] for the second term
        I=(I_0*math.exp((-q*V)/(n*k*T))*(math.exp((q*V)/(k*T))-1))+((9/8)*p_s*p_0*mu_p*(V**2/(L)**m))
        
        
    if V>V_t:
    
       I=(I_0*math.exp((-q*V)/(n*k*T))*(math.exp((q*V)/(k*T))-1))+((9/8)*p_s*p_0*mu_p*(V**2/(L)**m))
    
# [5] also discusses a third linear region with V>V_er
    #if V>V_er: 
        #I=v/R_int
       
    V_arr.append(V)
    I_arr.append(I)

print("Voltage array")    
print(V_arr)
print("")
print("Current Array")
print(I_arr)

# grid lines
plt.grid(True)
plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

# add the labels you wish for x and y axis
plt.xlabel('Voltage (V)')
plt.ylabel('Current Density (A/m^2)')


#Plot
plt.plot(V_arr,I_arr)
#plt.semilogy(V_arr,I_arr)

plt.show()




# References:
#[1] A. Traore, “High Power Diamond Schottky Diode,” Univ. Grenoble, 2014.
#[2] [1] P. E. H. Rhoderick, M. Sc, D. Ph, C. Eng, and F. P, “Metal-semiconductor contacts,” IEE Rev., vol. 129, no. 1, pp. 1–14, 1982.
#[3]V. W. L. Chin, J. W. V. Storey, and M. A. Green, “P-type PtSi Schottky-diode barrier height determined from I-V measurement,” 
#Solid State Electron., vol. 32, no. 6, pp. 475–478, 1989.
#[4] V. W. L. Chin, J. W. V. Storey, and M. A. Green, “Characteristics of p-type PtSi Schottky diodes under reverse bias,” 
#J. Appl. Phys., vol. 68, no. 8, pp. 4127–4132, 1990.
#[5]M. Brezeanu, “Diamond Schottky Barrier Diodes (Doctoral Thesis),” University of Cambridge, 2007.
