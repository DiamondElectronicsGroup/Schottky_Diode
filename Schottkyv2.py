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
#p_0 is the permitivity of free space (m^-3kg^-1s^4A^2)
p_0=8.85E-12
#pi
pi=math.pi

#Material Parameters
#METAL
#wf_m is the work function (eV)
wf_m=3.8

#SEMICONDUCTOR 
#The following parameters are for Boron-doped diamond
#Eg is the Bandgap(eV)
Eg=5.47
#xi is the difference between the Fermi level and the Valence band (eV)[5]
xi=0.3
#E_max is the maximum electric field in reverse bias for IFL [7] 0.8-1.4 MVcm^_1 = 1E8v/m
E_max=1E8
#wf_s is the work function(eV) [5]
wf_s=4.2
#X_s is the electron affinity(eV) (positive or negative??)[5]
X_s=1
#Na is the doping concentration m^-3 (between 10^20 and 10^24)[5]
Na=1E22
#p_s is the relative permitivity[5]
p_s=6.21
#mu_p is the hole mobility(m^2 V^-1 s^-1)
mu_p=0.2
#Internal reistance of the semiconductor (ohms)
R_int=4.85E-10
#R is the effective Richardson constant A m^-2 K^-2 given in [1,6]
R=90.0E4
#n is the ideality factor [6]
n=1.15
#n_i is the intrinsic carrier density (m^-3)
n_i=6E13
#tau s the lifetime of carriers in the depletion region(s)
tau=7E-9

#OTHER PARAMETERS
#Sizes
#A is the diode area m^2 
#(the size of the schottky contact pi*90nm^2 in m^2)
A=2.54E-14
#L is the schottky diode channel length in m
L= 20E-9

#T is the absolute temperature i.e. room temperature (K)
T=296.0 

#----------------------------------------------------
#CALCULATIONS

#SB_0 is the Schottky barrier at zero bias [2]
SB_0=Eg-(wf_m-X_s)

#V_bi is the built in voltage (eV)
V_bi=wf_s-wf_m

#m_e is the effetive mass of holes in the semiconductor
m_e=0.4*me

#Electron Tunneling Current density parameters[4]:
E_inf=(q*h/(4*pi))*math.sqrt(Na/(m_e*p_s*p_0))
E_0=E_inf/math.tanh(E_inf/(k*T))
 
#----------------------------------------------------------
#DETERMINATION OF THE CURRENT

#There are various transport mechanisms which determine the current density of the Schottky diode


J=0.0
V=0.0
V_arr=[]
J_arr=[]

SB=0.0

J_TED0=0.0
J_TED=0.0

J_MG=0.0

J_tun0=0.0
J_tun1=0.0
J_tun=0.0

W=0.0
J_r=0.0

J_IFL0=0.0
J_IFL1=0.0
J_IFL=0.0

J_TFE0=0.0
J_TFE1=0.0
J_TFE2=0.0
J_TFE=0.0

for i in range (-1001,1001):

    
    V=i/100.0
    SB=SB_0+((1/n)*V)
    

    # Current Densities derived in [4]  
    #Thermionic emission Diffusion Current Density [1,2,4,7]
    J_TED0=R*T**2*math.exp((-1*q*SB)/(k*T))
    J_TED=abs(J_TED0*(math.exp((q*V)/(k*T))-1))
    J_TED=0

    #Space Charge Limited Current Density[7]
    #J_MG=((9/8)*p_s*p_0*mu_p*(V**2/(L)**3))
    #J_MG=0  
        
    #Recombination Current Density[2,4]
    W=((2*p_s*p_0*(V_bi-V))/(q*Na))**0.5
    J_r=abs(((q*n_i*W)/((V_bi-V)*tau))*(math.exp((q*V)/(2*k*T))-1))
    #J_r=0    

    #Quantum Tunneling [4]
    
    J_tun0=(R*T*math.sqrt(pi*E_inf*q*(SB-V-xi)))/(k*(1/math.tanh(E_inf/k*T)))
    J_tun1=math.exp((-1*q*xi)/(k*T))*math.exp((-1*q*(SB-xi))/E_0)  
    J_tun=J_tun0*J_tun1*math.exp(q*V/E_0)
    J_tun=0
        

       
    #Current Densities derived in [6]:
    #Image-Force Lowering
    J_IFL0=(-1*q*SB)/(k*T)
    J_IFL1=(q/(k*T))*math.sqrt((q*E_max)/(4*pi*p_s*p_0))
    J_IFL=R*T**2*math.exp(J_IFL0)*math.exp(J_IFL1)
    J_IFL=0 
    
    #Thermionic-field Emission
    #J_TFE0=(R*T*math.sqrt(pi*E_inf))/k
    #J_TFE1=(q*(V-xi)+(q*SB)/(math.cosh((E_inf)/(k*T))**2))**0.5
    #J_TFE2=math.exp((-1*SB)/E_0)
    #J_TFE=math.abs(J_TFE0*J_TFE1*J_TFE2*math.exp(q*V*((1/(k*T))-(1/E_0))))
    #J_TFE=0
    
    #Total Current
    J=J_TED+J_MG+J_tun+J_r+J_IFL+J_TFE
    

    V_arr.append(V)
    J_arr.append(J)


# grid lines
plt.grid(True)
plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

# add the labels you wish for x and y axis
plt.xlabel('Voltage (V)')
plt.ylabel('Current Density (A/m^2)')


#Plot
plt.plot(V_arr,J_arr)
plt.semilogy(V_arr,J_arr)

plt.show()

#---------------------------------------------------------------------------------------------------------------------------------------
# References:
#[1] A. Traore, “High Power Diamond Schottky Diode,” Univ. Grenoble, 2014.
#[2] [1] P. E. H. Rhoderick, M. Sc, D. Ph, C. Eng, and F. P, “Metal-semiconductor contacts,” IEE Rev., vol. 129, no. 1, pp. 1–14, 1982.
#[3]V. W. L. Chin, J. W. V. Storey, and M. A. Green, “P-type PtSi Schottky-diode barrier height determined from I-V measurement,” 
#Solid State Electron., vol. 32, no. 6, pp. 475–478, 1989.
#[4] S. Amor, A. Ahaitouf, P. Salvestrini, J, and A. Ougazzaden, “Transport mechanisms in Schottky diodes realized on GaN,” 
#IOP Conf. Ser. Mater. Sci. Eng., vol. 186, no. 1, 2017.
#[5] L. Diederich, O. M. Küttel, P. Aebi, and L. Schlapbach, “Electron affinity and work function of differently oriented and doped diamond 
#surfaces determined by photoelectron spectroscopy,” Surface Science, vol. 418, no. 1. pp. 219–239, 1998.
#[6]T. Teraji, A. Fiori, N. Kiritani, S. Tanimoto, E. Gheeraert, and Y. Koide, “Mechanism of reverse current increase of vertical-type 
#diamond Schottky diodes,” J. Appl. Phys., vol. 122, no. 13, 2017.
#[7]I. Ben-Yaacov and U. K. Mishra, “Unipolar Space Charge Limited Transport,” in Univ. Of California Santa Barbara lecture notes, 2004, pp. 1–21.