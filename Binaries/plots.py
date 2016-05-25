## This code make scatter plots and histograms for the population of Binary neutron stars generated with BSE.
## Last modified 25th May.

#Inputs: "Binaries.dat", which is the bns.py output.
#Outputs: ".png" plots

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit

#Read datafiles

data = np.genfromtxt('binariaes.dat')

j = data[:,0]
i = data[:,1]
mass1 = data[:,2]
mass2 = data[:,3]
mns1 = data[:,4]
mns2 = data[:,5]
sntime1 = data[:,6]
sntime2 = data[:,7]
e0 = data[:,8]
a0 = data[:,9]
e = data[:,10]
a = data[:,11]
t = data[:,12]
vx1 = data[:,13]
vy1 = data[:,14]
vz1 = data[:,15]
vx2 = data[:,16]
vy2 = data[:,17]
vz2 = data[:,18]
z = data[:,19]

#Magnitude of center mass velocity after kick 1 and 2
v1 = (vx1**2 + vy1**2 + vz1**2)**(0.5)
v2 = (vx2**2 + vy2**2 + vz2**2)**(0.5)

############################General plots

### Neutron star mass vs Initial mass, dots with different colors, due to metallicity.

#####Conditions consider different intervals of metallicity
###Condition 1
ns1=[]
m1=[]
a1 = 4e-4
b1 = 7.8e-3

##Condition 2
ns2=[]
m2=[]
a2 = 7.8e-3
b2 = 1.5e-2

## Condition 3
ns3=[]
m3=[]
a3 = 1.5e-2
b3 = 2.2e-2 

##Condition 4
ns4=[]
m4=[]
a4 = 2.2e-2
b4 = 3e-2 

for i in range(len(t)):
     
        if (z[i] > a1 and z[i] < b1 ):
            ns1.append(mns1[i])
            m1.append(mass1[i])
        if (z[i] >= a2 and z[i] < b2 ):
            ns2.append(mns1[i])
            m2.append(mass1[i])
        if (z[i] >= a3 and z[i] < b3):
            ns3.append(mns1[i])
            m3.append(mass1[i])
        if (z[i] >= a4 and z[i] < b4):
            ns4.append(mns1[i])
            m4.append(mass1[i])
            
plt.figure(0)
line1,=plt.plot(m1,ns1,'b.',label='%.f <Z< %.f'%(a1,b1))
line2,=plt.plot(m2,ns2,'g.',label='%.f <Z< %.f'%(a2,b2))
line3,=plt.plot(m3,ns3,'r.',label='%.f <Z< %.f'%(a3,b3))
line4,=plt.plot(m4,ns4,'k.',label='%.f <Z< %.f'%(a4,b4))
plt.legend(handles=[line1,line2,line3,line4], loc=1)    
plt.xlabel('Initial mass (Msun)')
plt.ylabel('Neutron star mass (Msun)')
plt.savefig("mns1_vs_mass1.png")


##### Decay timescale vs Initial mass
plt.figure(2)
plt.plot(mass1,t,'b.')
plt.xlabel('Initial mass (Msun)')
plt.ylabel('Decay timescale (Myr)')
plt.title('Decay timescale  vs Mass')
plt.savefig("t_vs_mass1.png")

plt.figure(3)
plt.plot(mass2,t,'b.')
plt.xlabel('Initial mass (Msun)')
plt.ylabel('Decay timescale (Myr)')
plt.title('Decay timescale vs Mass')
plt.savefig("t_vs_mass2.png")

###### Supernova time vs Metallicity
plt.figure(4)
plt.plot(z,sntime1,'b.')
plt.xlabel('Metalicity')
plt.ylabel('supernova time (Myr)')
plt.title('Supernova time vs Metallicity ')
plt.savefig("sntime1_vs_z.png")

plt.figure(5)
plt.plot(z,sntime2,'b.')
plt.xlabel('Metalicity')
plt.ylabel('supernova time 2 (Myr)')
plt.title('Supernova time vs Metallicity ')
plt.savefig("sntime2_vs_z.png")


##### Neutron star mass vs Metallicity
plt.figure(6)
plt.plot(z,mns1,'b.')
plt.xlabel('Metalicity')
plt.ylabel('Neutron star mass')
#plt.title('Neutron star mass vs Metallicity ')
plt.savefig('mns1_vs_z.png')


plt.figure(7)
plt.plot(z,mns2,'b.')
plt.xlabel('Metalicity')
plt.ylabel('Neutron star mass')
#plt.title('Neutron star mass vs Metallicity')
plt.savefig('mns2_vs_z.png')



######## Decay timescale vs semimajor axis
plt.figure(8)
plt.plot(a,t,'b.')
plt.xlabel('Semimajor axis')
plt.xscale('log')
plt.ylabel('Decay timescale')
plt.title('Decay timescale vs semimajor axis')
plt.savefig('t_vs_a.png')


####### Decay timescales vs Eccentricity
plt.figure(9)
plt.plot(e,t,'b.')
plt.xlabel('Eccentricity')
plt.ylabel('Decay timescale')
plt.title('Decay timescale vs Eccentricity')
plt.savefig('t_vs_e.png')

###### CM velocity 2 vs CM velocity 1
plt.figure(10)
plt.plot(v1,v2,'b.')
plt.xlabel('V_cm kick1')
plt.ylabel('V_cm kick2')
plt.savefig('Vcm2_vs_Vcm1.png')


######## Supernova time vs Total mass

####### Conditions to different metallicity intervals
###Condition 1
sn1=[]
M1=[]
a1 = 4e-4
b1 = 7.8e-3

##Condition 2
sn2=[]
M2=[]
a2 = 7.8e-3
b2 = 1.5e-2

## Condition 3
sn3=[]
M3=[]
a3 = 1.5e-2
b3 = 2.2e-2 

##Condition 4
sn4=[]
M4=[]
a4 = 2.2e-2
b4 = 3e-2 

i = 0


for i in range(len(mass1)):
    M = mass1[i] + mass2[i]

    if mass1[i]>mass2[i]:
        
        if (z[i] > a1 and z[i] < b1 ):
            sn1.append(sntime2[i])
            M1.append(M)
        if (z[i] >= a2 and z[i] < b2 ):
            sn2.append(sntime2[i])
            M2.append(M)
        if (z[i] >= a3 and z[i] < b3 ):
            sn3.append(sntime2[i])
            M3.append(M)
        if (z[i] >= a4 and z[i] < b4 ):
            sn4.append(sntime2[i])
            M4.append(M)
            
plt.figure(11)
line1,=plt.plot(M1,sn1,'b.',label='%.f <Z< %.f'%(a1,b1))
line2,=plt.plot(M2,sn2,'g.',label='%.f <Z< %.f'%(a2,b2))
line3,=plt.plot(M3,sn3,'r.',label='%.f <Z< %.f'%(a3,b3))
line4,=plt.plot(M4,sn4,'k.',label='%.f <Z< %.f'%(a4,b4))
plt.legend(handles=[line1,line2,line3,line4], loc=1)
plt.xlabel('Total mass (MSun)')
plt.ylabel('Supernova time (Myr)')
plt.savefig('supernovae.png')



#####Histograms 

plt.figure(13)
bins = np.arange(0, 1000, 30)
plt.hist(a, bins)
plt.xlabel('semimajor axis')
plt.ylabel('frequency')
plt.title('Semimajor axis')
plt.savefig('Histogram_semimajor.png')

plt.figure(14)
plt.hist(e) 
plt.xlabel('eccentricity')
plt.ylabel('frequency')
plt.title('Histogram eccentricity')
plt.savefig('Histogram_eccentricity.png')

plt.figure(15)
plt.hist(v2)
plt.xlabel('Second kick velocity (Km/s)')
plt.ylabel('frequency')
plt.title('Center of mass velocity')
plt.savefig('Histogram_vcm2.png')

plt.figure(16)
plt.hist(v1)
plt.xlabel('Kick velocity (Km/s)')
plt.ylabel('frequency')
plt.title('Kick velocity')
plt.savefig('Histogram_vcm1.png')

plt.figure(17)
plt.hist(t)
plt.xlabel('Decay timsecale (Myr)')
plt.ylabel('Frequency')
plt.title('Decay timescale')
plt.savefig('Histogram_decay_timescale.png')



###### Decay timescale fit to a power law
plt.figure(18)
        
tmin=0.0
tmax=3000.0
dt=100.0
bins = np.arange(tmin,tmax,dt)
t=t[(t<tmax)*(t>=tmin)]
freq,bins = np.histogram(t,bins)
#xdata = (bins[:-1]+bins[1:])/2
xdata = bins[:-1]
#ydata = freq[:]

i=2  
ydata=np.zeros(len(freq))
for i in xrange(len(freq)):
    ydata[i]=sum(freq[:i])
    

def f(t,t0,k,N0):
    return N0*(t/t0)**(1-k)

popt,pcov =curve_fit(f,xdata,ydata,[100,0.2,100])
plt.plot(xdata,ydata,'k.')
xfit=np.linspace(xdata[0],xdata[-1],1000)
plt.plot(xfit,f(xfit,popt[0],popt[1],popt[2]))
plt.xlabel('Decay timescale (Myr)')
plt.ylabel('frequency')
plt.title('Decay timescale (Myr) ; N0 = %0.2f, t0=%.3f, k=%.3f'%(popt[2],popt[0],popt[1]))
plt.savefig('t_decay_fit_log.png')



##### Final vs Initial semimajor axis. Color are eccentricities

###Conditions for different eccentricity intervals
#Condition 1
a1=[]
a01=[]
e1 = 0.0
e2 = 0.2

#Condition2
a2=[]
a02=[]
e3 = 0.2
e4 = 0.4

#Condition 3
a3=[]
a03=[]
e5 = 0.4
e6 = 0.6

#Condition 4
a4=[]
a04=[]
e7 = 0.6
e8 = 0.8

#Condition 5
a5=[]
a05=[]
e9 = 0.8
e10 = 1.0

for k in range(len(t)):
    if (e[k] > e1 and e[k] < e2 ):
        a1.append(a[k])
        a01.append(a0[k])
    if (e[k] >= e3 and e[k] < e4 ):
        a2.append(a[k])
        a02.append(a0[k])
    if (e[k] >= e5 and e[k] < e6 ):
        a3.append(a[k])
        a03.append(a0[k])
    if (e[k] >= e7 and e[k] < e8 ):
        a4.append(a[k])
        a04.append(a0[k])
    if (e[k] >= e9 and e[k] < e10 ):
        a5.append(a[k])
        a05.append(a0[k])
    

plt.figure(19)
line1,=plt.plot(a01,a1,'b.',label='0.0<e<0.2')
line2,=plt.plot(a02,a2,'r.',label='0.2<e<0.4')
line3,=plt.plot(a03,a3,'g.',label='0.4<e<0.6')
line4,=plt.plot(a04,a4,'k.',label='0.6<e<0.8')
line5,=plt.plot(a05,a5,'y.',label='0.8<e<1.0')
plt.legend(handles=[line1,line2,line3,line4,line5], loc=1)
plt.xlabel('Initial semimajor axis (Rsun)')
plt.ylabel('Final semimajor axis (Rsun)')
plt.yscale('log')
plt.savefig('separations_e.png')


####### Semimajor axis vs eccentricity
plt.figure(20)
plt.plot(e,a,'k.')
plt.ylim(0,1000)
plt.savefig('dependence.png')


####### Final vs Initial eccentricity
plt.figure(20)
plt.plot(e0,e,'b.')
plt.xlabel('Initial eccentricity')
plt.ylabel('Final eccentricity')
#plt.title('Initial vs final eccentricity')
plt.savefig('eccentricities.png')


#######Histogram Center of mass velocicities
plt.figure(21)
bins=np.arange(0,500,20)
plt.hist(v2,bins)
plt.hist(v1,bins)
plt.xlabel('Center of mass Velocity (Km/s)')
plt.ylabel('frequency')
plt.savefig('Histograms_kick.png')


####### Histogram initial eccentricities
plt.figure(22)
plt.hist(e0)
plt.title('Histogram initial eccentricity')
plt.xlabel('Initial eccentricity')
plt.ylabel('frequency')
plt.savefig('Histogram_e0.png') 

####### Eccentricity-period distribution

grav=6.67e-11
Msun=2e30
Rsun=6.96e8
convert = 1.15741e-5

period=[]
for i in range(len(e)):
    M=(mns1[i]+mns2[i])*Msun
    T = ( ((4*np.pi*np.pi) / (grav*M)) *(a[i]*Rsun)**3 )**0.5
    T = T*convert
    print mns1[i], mns2[i], a[i], T
    period.append(T)

plt.figure(0)
plt.plot(period,e, 'b.')
plt.xlim(0.1,100)
plt.xscale('log')
plt.ylabel('Eccentricity')
plt.xlabel('Period (days)')
plt.savefig('period_distribution1.png')

plt.figure(23)
bins = np.arange(0,150,10)
plt.hist(period,bins)
plt.xlabel('Period (days)')
plt.savefig('period_hist.png')

