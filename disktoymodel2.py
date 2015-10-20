# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:57:00 2015

@author: jessicaluna
"""

from pylab import*
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import csv
import glob
import time
import scipy.integrate as integrate

starttime=time.time()


###PART 1: CALCULATING FLUX DENISTY AT SPECITFIED ORBITAL RADIUS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
def fluxdensity(star_radius,Temp,orbital_radius,freq):
   ##DEFINING CONSTANTS
    h=6.626*10**(-27) #erg s^-1
    c=3.0*10**(10) #cm s^-1
    k_b= 1.38*10**(-16)  #
    
    B_v=(2.0*h*pow(freq,3) /(pow(c,2)))*(1/(exp((h*freq)/(k_b*Temp))-1.0))
    F_v=B_v *(pi)*(10**(23))
    
    Flux_absorbed_lis=((star_radius/orbital_radius)**2) *F_v
    return (Flux_absorbed_lis) 

c_speed=3.0*10**(10) #cm s^-1
wave_lis,Q_abs_01=loadtxt('01micronmodel.txt',unpack=True,usecols=[0,1])
Q_abs_1=loadtxt('1micronmodel.txt',unpack=True,usecols=[1])
Q_abs_10=loadtxt('10micronmodel.txt',unpack=True,usecols=[1])

wave=wave_lis*10**(-4)
#print(wave)

freq_lis=c_speed/wave

#print(freq_lis)

R_sun= 6.96*10**(10) #cm
R_star_fom= 1.842*R_sun
Temp_fom= 8590 #K 
AU_CM=1.496*10**(13) ##cm in 1AU
radius_10AU=10.0*AU_CM
radius_130AU=130.0*AU_CM


Flux_absorbed_lis_10AU = fluxdensity(R_star_fom,Temp_fom,radius_10AU ,freq_lis)
Flux_absorbed_lis_130AU = fluxdensity(R_star_fom,Temp_fom, radius_130AU,freq_lis)



##part 2 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
def integrateflux(a,freq,flux,Q_abs):
    flux=flux*10**(-23)
    integrand= pi*(a**2)*flux*Q_abs
    Pin=integrate.trapz(integrand,freq)
    return Pin
    
a_01=0.1*10**(-4)
a_1=1.0*10**(-4)
a_10=10.0*10**(-4)

##calulating total energy absorbed by dust grain in erg/s!
Pin_01_10AU=integrateflux(a_01,freq_lis, Flux_absorbed_lis_10AU ,Q_abs_01)
Pin_01_130AU=integrateflux(a_01,freq_lis, Flux_absorbed_lis_130AU ,Q_abs_01)

Pin_1_10AU=integrateflux(a_1,freq_lis, Flux_absorbed_lis_10AU ,Q_abs_1)
Pin_1_130AU=integrateflux(a_1,freq_lis, Flux_absorbed_lis_130AU ,Q_abs_1)

Pin_10_10AU=integrateflux(a_10,freq_lis, Flux_absorbed_lis_10AU ,Q_abs_10)
Pin_10_130AU=integrateflux(a_10,freq_lis, Flux_absorbed_lis_130AU ,Q_abs_10)
print('Pin 0.1u at 10aU',Pin_01_10AU)#*10**(-7))
print('Pin 0.1u at 130aU',Pin_01_130AU)#*10**(-7))
print('Pin 1u at 10aU',Pin_1_10AU)#*10**(-7))
print('Pin 1u at 130aU',Pin_1_130AU)#*10**(-7))
print('Pin 10u at 10aU',Pin_10_10AU)#*10**(-7))
print('Pin 10u at 130aU',Pin_10_130AU)#*10**(-7))


##part 4: calculating equilibrium temperature and emission spectrum;;;;;;;;;;;;;;;;;;;;;;;
#calculating Pout of a dust grain for an arbitrary temp
def pout(Temp,Qabs,a ,freq):
   ##DEFINING CONSTANTS
    h=6.626*10**(-27) #erg s^-1
    c=3.0*10**(10) #cm s^-1
    k_b= 1.38*10**(-16)  #
    
    B_v=(2.0*h*pow(freq,3) /(pow(c,2)))*(1/(exp((h*freq)/(k_b*Temp))-1.0))
    F_v=B_v *(pi)
    integrand= 4*pi*(a**2)* F_v*Qabs
    Pout=integrate.trapz(integrand,freq)
    return (Pout)

#finding the temp of dust grain by setting it equal to pin

max_error=1.0*10**(-4)

#check starting temp-- blackbody ;;;;;;;;;;;;;;;;;;;;
def equil_blackbody(Pin,ass):
    sigma=5.67*10**(-5) 
    Tin=(Pin/(4*pi*sigma*(ass**2)))**(0.25)   
    return Tin,Tout
    
def bisection_temp(T_start,T_end,tol,Pin,Q_abs,a,freq_lis):
    Tc=(T_end+T_start)/2.0
    n=1
    n_max=100
    Pc=pout(Tc,Q_abs,a,freq_lis)
    
    while (abs(Pc-Pin)/Pin > tol) & (n<n_max) :       
        if Pc>Pin:
            T_end=Tc
        else :
            T_start = Tc
        Tc = (T_end+T_start)/2.0
        #print(n,Tc,Pin,Pc)
        n=n+1
        Pc=pout(Tc,Q_abs,a,freq_lis)
 
    return Tc

Tstart=1.0
Tend= 1.0*10**(5)
T_equil_01_10AU= bisection_temp(Tstart,Tend,max_error,Pin_01_10AU,Q_abs_01,a_01,freq_lis) 
T_equil_01_130AU= bisection_temp(Tstart,Tend,max_error,Pin_01_130AU,Q_abs_01,a_01,freq_lis)  
T_equil_1_10AU= bisection_temp(Tstart,Tend,max_error,Pin_1_10AU,Q_abs_1,a_1,freq_lis) 
T_equil_1_130AU= bisection_temp(Tstart,Tend,max_error,Pin_1_130AU,Q_abs_1,a_1,freq_lis)
T_equil_10_10AU= bisection_temp(Tstart,Tend,max_error,Pin_10_10AU,Q_abs_10,a_10,freq_lis)  
T_equil_10_130AU= bisection_temp(Tstart,Tend,max_error,Pin_10_130AU,Q_abs_10,a_10,freq_lis)   


  
print('T_equil 0.1u 10AU=',T_equil_01_10AU,pout(T_equil_01_10AU,Q_abs_01,a_01,freq_lis))
print('T_equil 0.1u 130=',T_equil_01_130AU)
print('T_equil 1u 10AU=',T_equil_1_10AU)
print('T_equil 1u 130=',T_equil_1_130AU)
print('T_equil 10u 10AU=',T_equil_10_10AU)
print('T_equil 10u 130=',T_equil_10_130AU)

   
def flux(Temp ,freq,a,Q_abs):
   ##DEFINING CONSTANTS
    h=6.626*10**(-27) #erg s^-1
    c=3.0*10**(10) #cm s^-1
    k_b= 1.38*10**(-16)  #

    B_v=(2.0*h*pow(freq,3) /(pow(c,2)))*(1/(exp((h*freq)/(k_b*Temp))-1.0))
    F_v=B_v *(pi)*10**(23)*4.0*pi*(a**2)*Q_abs
    return(F_v)

equilFlux_01_10AU=flux(T_equil_01_10AU,freq_lis,a_01, Q_abs_01)
equilFlux_01_130AU=flux(T_equil_01_130AU,freq_lis,a_01,Q_abs_01)    
equilFlux_1_10AU=flux(T_equil_1_10AU,freq_lis,a_1, Q_abs_1) 
equilFlux_1_130AU=flux(T_equil_1_130AU,freq_lis,a_1, Q_abs_1) 
equilFlux_10_10AU=flux(T_equil_1_10AU,freq_lis,a_10, Q_abs_10) 
equilFlux_10_130AU=flux(T_equil_10_130AU,freq_lis,a_10, Q_abs_10) 


###part 4 finding the number of dust grains and dust mass
flux_130AU=10.0  ##peaks at 10 AU in janskys
flux_10AU=1.0  ##peak at 10 AU in janskys

def dustmass(equilflux,a,flux_scale):
    pc_cm=3.086*10**(18) #cm in pc
    dist=7.61*pc_cm
    peakflux=max(equilflux)*(1/dist)**2  
    
    density=2.0 # g/cm^3

    Vol=(4.0/3.0)*pi*(a**3)
    n=flux_scale/peakflux
    mass=n*density*Vol
    return (n,mass)
    

n_01_10AU,mass_01_10AU=dustmass(equilFlux_01_10AU,a_01,flux_10AU)
n_01_130AU,mass_01_130AU=dustmass(equilFlux_01_130AU,a_01,flux_130AU)
n_1_10AU,mass_1_10AU=dustmass(equilFlux_1_10AU,a_1,flux_10AU)
n_1_130AU,mass_1_130AU=dustmass(equilFlux_1_130AU,a_1,flux_130AU)
n_10_10AU,mass_10_10AU=dustmass(equilFlux_10_10AU,a_10,flux_10AU)
n_10_130AU,mass_10_130AU=dustmass(equilFlux_10_130AU,a_10,flux_130AU)


print('n          ','             Mass')
print(n_01_10AU,mass_01_10AU,mass_01_10AU/n_01_10AU)
print(n_01_130AU,mass_01_130AU)
print(n_1_10AU,mass_1_10AU)
print(n_1_130AU,mass_1_130AU)
print(n_10_10AU,mass_10_10AU)
print(n_10_130AU,mass_10_130AU)
    

##part 5 calculating the radiation pressure and Poynting-Robertson drag and timescale
def radpressure_PRdrag(Pin,orbital_radius,mass,n):
    md=mass/n
    c=3.0*10**(10) #cm s^-1
    G=6.67*10**(-8) ##cm3 g-1 s-2
    M_sun=1.989*10**(33) #g
    Ms=1.92*M_sun
    F_radiation=Pin/c
    v=sqrt(G*Ms/orbital_radius)
    F_pr=v*Pin/(c**2)
    F_grav=G*Ms*md/(orbital_radius**2)

    beta=F_radiation/F_grav
    timescale=10**(-23)*400.0*(orbital_radius**2)/beta
    
    
    return (F_radiation,F_pr,timescale)
    

 
F_rad_01_10AU,F_PR_01_10AU,t_01_10AU=radpressure_PRdrag(Pin_01_10AU,radius_10AU,mass_01_10AU,n_01_10AU) 
F_rad_01_130AU,F_PR_01_130AU,t_01_130AU=radpressure_PRdrag(Pin_01_130AU,radius_130AU,mass_01_130AU,n_01_130AU) 
F_rad_1_10AU,F_PR_1_10AU,t_1_10AU=radpressure_PRdrag(Pin_1_10AU,radius_10AU,mass_1_10AU,n_1_10AU) 
F_rad_1_130AU,F_PR_1_130AU,t_1_130AU=radpressure_PRdrag(Pin_1_130AU,radius_130AU,mass_1_130AU,n_1_130AU) 
F_rad_10_10AU,F_PR_10_10AU,t_10_10AU=radpressure_PRdrag(Pin_10_10AU,radius_10AU,mass_10_10AU,n_10_10AU) 
F_rad_10_130AU,F_PR_10_130AU,t_10_130AU=radpressure_PRdrag(Pin_10_130AU,radius_130AU,mass_10_130AU,n_10_130AU)   

sec_yr=3.15*10**(7)
print('F_rad, F_pr','timescale')
print(F_rad_01_10AU,F_PR_01_10AU,t_01_10AU/(sec_yr) )
print(F_rad_1_10AU,F_PR_1_10AU,t_1_10AU/sec_yr) 
print(F_rad_10_10AU,F_PR_10_10AU,t_10_10AU/sec_yr)
print(F_rad_01_130AU,F_PR_01_130AU,t_01_130AU/sec_yr)
print(F_rad_1_130AU,F_PR_1_130AU,t_1_130AU/sec_yr)
print(F_rad_10_130AU,F_PR_10_130AU,t_10_130AU/sec_yr)         


figure(1)
plt.ylabel('Absorbed log(Flux) Janskys')
plt.xlabel('log (wavelength) um')
plt.title('Flux Density')
plt.loglog(wave_lis[0:190],equilFlux_01_10AU[0:190],label='0.1 um 10AU')
plt.loglog(wave_lis[0:190],equilFlux_01_130AU[0:190],label='0.1 um 130AU')
plt.loglog(wave_lis[0:190],equilFlux_1_10AU[0:190],label='1 um 10AU')
plt.loglog(wave_lis[0:190],equilFlux_1_130AU[0:190],label='1 um 130AU')
plt.loglog(wave_lis[0:190],equilFlux_10_10AU[0:190],label='10 um 10AU')
plt.loglog(wave_lis[0:190],equilFlux_10_130AU[0:190],label='10 um 130AU')
plt.axis([min(wave_lis),10**4,10**(-4),10**(11)])

plt.legend( loc=2)

figure(2)
plt.ylabel('Absorbed log(Flux) Janskys')
plt.xlabel('log (wavelength) um')
plt.title('Flux Density')
plt.loglog(wave_lis[0:185]*10**4,Flux_absorbed_lis_10AU[0:185],c='r',label='10 AU')
plt.loglog(wave_lis[0:185]*10**4,Flux_absorbed_lis_130AU[0:185],label='130 AU')
plt.legend( loc=4)

show()    
    







print('Time to run',time.time()-starttime)
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:57:00 2015

@author: jessicaluna
"""

from pylab import*
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import csv
import glob
import time
import scipy.integrate as integrate

starttime=time.time()


###PART 1: CALCULATING FLUX DENISTY AT SPECITFIED ORBITAL RADIUS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
def fluxdensity(star_radius,Temp,orbital_radius,freq):
   ##DEFINING CONSTANTS
    h=6.626*10**(-27) #erg s^-1
    c=3.0*10**(10) #cm s^-1
    k_b= 1.38*10**(-16)  #
    
    B_v=(2.0*h*pow(freq,3) /(pow(c,2)))*(1/(exp((h*freq)/(k_b*Temp))-1.0))
    F_v=B_v *(pi)*(10**(23))
    
    Flux_absorbed_lis=((star_radius/orbital_radius)**2) *F_v
    return (Flux_absorbed_lis) 

c_speed=3.0*10**(10) #cm s^-1
wave_lis,Q_abs_01=loadtxt('01micronmodel.txt',unpack=True,usecols=[0,1])
Q_abs_1=loadtxt('1micronmodel.txt',unpack=True,usecols=[1])
Q_abs_10=loadtxt('10micronmodel.txt',unpack=True,usecols=[1])

wave=wave_lis*10**(-4)
#print(wave)

freq_lis=c_speed/wave

#print(freq_lis)

R_sun= 6.96*10**(10) #cm
R_star_fom= 1.842*R_sun
Temp_fom= 8590 #K 
AU_CM=1.496*10**(13) ##cm in 1AU
radius_10AU=10.0*AU_CM
radius_130AU=130.0*AU_CM


Flux_absorbed_lis_10AU = fluxdensity(R_star_fom,Temp_fom,radius_10AU ,freq_lis)
Flux_absorbed_lis_130AU = fluxdensity(R_star_fom,Temp_fom, radius_130AU,freq_lis)



##part 2 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
def integrateflux(a,freq,flux,Q_abs):
    flux=flux*10**(-23)
    integrand= pi*(a**2)*flux*Q_abs
    Pin=integrate.trapz(integrand,freq)
    return Pin
    
a_01=0.1*10**(-4)
a_1=1.0*10**(-4)
a_10=10.0*10**(-4)

##calulating total energy absorbed by dust grain in erg/s!
Pin_01_10AU=integrateflux(a_01,freq_lis, Flux_absorbed_lis_10AU ,Q_abs_01)
Pin_01_130AU=integrateflux(a_01,freq_lis, Flux_absorbed_lis_130AU ,Q_abs_01)

Pin_1_10AU=integrateflux(a_1,freq_lis, Flux_absorbed_lis_10AU ,Q_abs_1)
Pin_1_130AU=integrateflux(a_1,freq_lis, Flux_absorbed_lis_130AU ,Q_abs_1)

Pin_10_10AU=integrateflux(a_10,freq_lis, Flux_absorbed_lis_10AU ,Q_abs_10)
Pin_10_130AU=integrateflux(a_10,freq_lis, Flux_absorbed_lis_130AU ,Q_abs_10)
print('Pin 0.1u at 10aU',Pin_01_10AU)#*10**(-7))
print('Pin 0.1u at 130aU',Pin_01_130AU)#*10**(-7))
print('Pin 1u at 10aU',Pin_1_10AU)#*10**(-7))
print('Pin 1u at 130aU',Pin_1_130AU)#*10**(-7))
print('Pin 10u at 10aU',Pin_10_10AU)#*10**(-7))
print('Pin 10u at 130aU',Pin_10_130AU)#*10**(-7))


##part 4: calculating equilibrium temperature and emission spectrum;;;;;;;;;;;;;;;;;;;;;;;
#calculating Pout of a dust grain for an arbitrary temp
def pout(Temp,Qabs,a ,freq):
   ##DEFINING CONSTANTS
    h=6.626*10**(-27) #erg s^-1
    c=3.0*10**(10) #cm s^-1
    k_b= 1.38*10**(-16)  #
    
    B_v=(2.0*h*pow(freq,3) /(pow(c,2)))*(1/(exp((h*freq)/(k_b*Temp))-1.0))
    F_v=B_v *(pi)
    integrand= 4*pi*(a**2)* F_v*Qabs
    Pout=integrate.trapz(integrand,freq)
    return (Pout)

#finding the temp of dust grain by setting it equal to pin

max_error=1.0*10**(-4)

#check starting temp-- blackbody ;;;;;;;;;;;;;;;;;;;;
def equil_blackbody(Pin,ass):
    sigma=5.67*10**(-5) 
    Tin=(Pin/(4*pi*sigma*(ass**2)))**(0.25)   
    return Tin,Tout
    
def bisection_temp(T_start,T_end,tol,Pin,Q_abs,a,freq_lis):
    Tc=(T_end+T_start)/2.0
    n=1
    n_max=100
    Pc=pout(Tc,Q_abs,a,freq_lis)
    
    while (abs(Pc-Pin)/Pin > tol) & (n<n_max) :       
        if Pc>Pin:
            T_end=Tc
        else :
            T_start = Tc
        Tc = (T_end+T_start)/2.0
        #print(n,Tc,Pin,Pc)
        n=n+1
        Pc=pout(Tc,Q_abs,a,freq_lis)
 
    return Tc

Tstart=1.0
Tend= 1.0*10**(5)
T_equil_01_10AU= bisection_temp(Tstart,Tend,max_error,Pin_01_10AU,Q_abs_01,a_01,freq_lis) 
T_equil_01_130AU= bisection_temp(Tstart,Tend,max_error,Pin_01_130AU,Q_abs_01,a_01,freq_lis)  
T_equil_1_10AU= bisection_temp(Tstart,Tend,max_error,Pin_1_10AU,Q_abs_1,a_1,freq_lis) 
T_equil_1_130AU= bisection_temp(Tstart,Tend,max_error,Pin_1_130AU,Q_abs_1,a_1,freq_lis)
T_equil_10_10AU= bisection_temp(Tstart,Tend,max_error,Pin_10_10AU,Q_abs_10,a_10,freq_lis)  
T_equil_10_130AU= bisection_temp(Tstart,Tend,max_error,Pin_10_130AU,Q_abs_10,a_10,freq_lis)   


  
print('T_equil 0.1u 10AU=',T_equil_01_10AU,pout(T_equil_01_10AU,Q_abs_01,a_01,freq_lis))
print('T_equil 0.1u 130=',T_equil_01_130AU)
print('T_equil 1u 10AU=',T_equil_1_10AU)
print('T_equil 1u 130=',T_equil_1_130AU)
print('T_equil 10u 10AU=',T_equil_10_10AU)
print('T_equil 10u 130=',T_equil_10_130AU)

   
def flux(Temp ,freq,a,Q_abs):
   ##DEFINING CONSTANTS
    h=6.626*10**(-27) #erg s^-1
    c=3.0*10**(10) #cm s^-1
    k_b= 1.38*10**(-16)  #

    B_v=(2.0*h*pow(freq,3) /(pow(c,2)))*(1/(exp((h*freq)/(k_b*Temp))-1.0))
    F_v=B_v *(pi)*10**(23)*4.0*pi*(a**2)*Q_abs
    return(F_v)

equilFlux_01_10AU=flux(T_equil_01_10AU,freq_lis,a_01, Q_abs_01)
equilFlux_01_130AU=flux(T_equil_01_130AU,freq_lis,a_01,Q_abs_01)    
equilFlux_1_10AU=flux(T_equil_1_10AU,freq_lis,a_1, Q_abs_1) 
equilFlux_1_130AU=flux(T_equil_1_130AU,freq_lis,a_1, Q_abs_1) 
equilFlux_10_10AU=flux(T_equil_1_10AU,freq_lis,a_10, Q_abs_10) 
equilFlux_10_130AU=flux(T_equil_10_130AU,freq_lis,a_10, Q_abs_10) 


###part 4 finding the number of dust grains and dust mass
flux_130AU=10.0  ##peaks at 10 AU in janskys
flux_10AU=1.0  ##peak at 10 AU in janskys

def dustmass(equilflux,a,flux_scale):
    pc_cm=3.086*10**(18) #cm in pc
    dist=7.61*pc_cm
    peakflux=max(equilflux)*(1/dist)**2  
    
    density=2.0 # g/cm^3

    Vol=(4.0/3.0)*pi*(a**3)
    n=flux_scale/peakflux
    mass=n*density*Vol
    return (n,mass)
    

n_01_10AU,mass_01_10AU=dustmass(equilFlux_01_10AU,a_01,flux_10AU)
n_01_130AU,mass_01_130AU=dustmass(equilFlux_01_130AU,a_01,flux_130AU)
n_1_10AU,mass_1_10AU=dustmass(equilFlux_1_10AU,a_1,flux_10AU)
n_1_130AU,mass_1_130AU=dustmass(equilFlux_1_130AU,a_1,flux_130AU)
n_10_10AU,mass_10_10AU=dustmass(equilFlux_10_10AU,a_10,flux_10AU)
n_10_130AU,mass_10_130AU=dustmass(equilFlux_10_130AU,a_10,flux_130AU)


print('n          ','             Mass')
print(n_01_10AU,mass_01_10AU,mass_01_10AU/n_01_10AU)
print(n_01_130AU,mass_01_130AU)
print(n_1_10AU,mass_1_10AU)
print(n_1_130AU,mass_1_130AU)
print(n_10_10AU,mass_10_10AU)
print(n_10_130AU,mass_10_130AU)
    

##part 5 calculating the radiation pressure and Poynting-Robertson drag and timescale
def radpressure_PRdrag(Pin,orbital_radius,mass,n):
    md=mass/n
    c=3.0*10**(10) #cm s^-1
    G=6.67*10**(-8) ##cm3 g-1 s-2
    M_sun=1.989*10**(33) #g
    Ms=1.92*M_sun
    F_radiation=Pin/c
    v=sqrt(G*Ms/orbital_radius)
    F_pr=v*Pin/(c**2)
    F_grav=G*Ms*md/(orbital_radius**2)

    beta=F_radiation/F_grav
    timescale=10**(-23)*400.0*(orbital_radius**2)/beta
    
    
    return (F_radiation,F_pr,timescale)
    

 
F_rad_01_10AU,F_PR_01_10AU,t_01_10AU=radpressure_PRdrag(Pin_01_10AU,radius_10AU,mass_01_10AU,n_01_10AU) 
F_rad_01_130AU,F_PR_01_130AU,t_01_130AU=radpressure_PRdrag(Pin_01_130AU,radius_130AU,mass_01_130AU,n_01_130AU) 
F_rad_1_10AU,F_PR_1_10AU,t_1_10AU=radpressure_PRdrag(Pin_1_10AU,radius_10AU,mass_1_10AU,n_1_10AU) 
F_rad_1_130AU,F_PR_1_130AU,t_1_130AU=radpressure_PRdrag(Pin_1_130AU,radius_130AU,mass_1_130AU,n_1_130AU) 
F_rad_10_10AU,F_PR_10_10AU,t_10_10AU=radpressure_PRdrag(Pin_10_10AU,radius_10AU,mass_10_10AU,n_10_10AU) 
F_rad_10_130AU,F_PR_10_130AU,t_10_130AU=radpressure_PRdrag(Pin_10_130AU,radius_130AU,mass_10_130AU,n_10_130AU)   

sec_yr=3.15*10**(7)
print('F_rad, F_pr','timescale')
print(F_rad_01_10AU,F_PR_01_10AU,t_01_10AU/(sec_yr) )
print(F_rad_1_10AU,F_PR_1_10AU,t_1_10AU/sec_yr) 
print(F_rad_10_10AU,F_PR_10_10AU,t_10_10AU/sec_yr)
print(F_rad_01_130AU,F_PR_01_130AU,t_01_130AU/sec_yr)
print(F_rad_1_130AU,F_PR_1_130AU,t_1_130AU/sec_yr)
print(F_rad_10_130AU,F_PR_10_130AU,t_10_130AU/sec_yr)         


figure(1)
plt.ylabel('Absorbed log(Flux) Janskys')
plt.xlabel('log (wavelength) um')
plt.title('Flux Density')
plt.loglog(wave_lis[0:190],equilFlux_01_10AU[0:190],label='0.1 um 10AU')
plt.loglog(wave_lis[0:190],equilFlux_01_130AU[0:190],label='0.1 um 130AU')
plt.loglog(wave_lis[0:190],equilFlux_1_10AU[0:190],label='1 um 10AU')
plt.loglog(wave_lis[0:190],equilFlux_1_130AU[0:190],label='1 um 130AU')
plt.loglog(wave_lis[0:190],equilFlux_10_10AU[0:190],label='10 um 10AU')
plt.loglog(wave_lis[0:190],equilFlux_10_130AU[0:190],label='10 um 130AU')
plt.axis([min(wave_lis),10**4,10**(-4),10**(11)])

plt.legend( loc=2)

figure(2)
plt.ylabel('Absorbed log(Flux) Janskys')
plt.xlabel('log (wavelength) um')
plt.title('Flux Density')
plt.loglog(wave_lis[0:185]*10**4,Flux_absorbed_lis_10AU[0:185],c='r',label='10 AU')
plt.loglog(wave_lis[0:185]*10**4,Flux_absorbed_lis_130AU[0:185],label='130 AU')
plt.legend( loc=4)

show()    
    







print('Time to run',time.time()-starttime)
