########### This code integrate the trajectory of binary systems under the Milky Way's Potential.
#Modified 11th May


from matplotlib.pyplot import *
from numpy import *
from scipy.integrate import odeint
from scipy.misc import derivative
from scipy import special
import numpy as np
from time import time as ptime
import os,commands
import pylab as p
import matplotlib.pyplot as plt
import random


##################Clean the directory

commands.getoutput(' rm orbits/orbits_*')
commands.getoutput(' rm orbits/orbits_general.dat')
commands.getoutput(' rm orbits/distrobution.dat')
commands.getoutput('rm orbits/times_*')
commands.getoutput('rm orbits/distribution.dat')


###################### Potential definition and differential ecuation systems definition
#Lenght units : 3.5 kpc
#Mass units : 5.5e10 Msun
#Time units : 

G = 1.0             #Gravitational constant
h = 1.0             #Lenght scale disk
M_disk = 1.0        #Mass of galactic disk
M_bulge = 0.21      #Bulge mass
M_halo = 17.85      #Halo mass
b = 0.2             #Vertical scale-lenght disk  
a = 0.34            #radial scale bulge
r_h = 9.96          #radial scale halo
r200 = 73.71        #Virial radius
rho_c = 1.12e-7     #Critical density
rho1 = 200*rho_c*(r200/r_h)*(1+r200/r_h)**2  
rho = M_bulge/2*np.pi*pow(a,3)

#### Definition of the  system of ordinary differential equiations. The solution are positions and velocities as function of time
def EOS(y,t):
    
    r=y[0:3]
    v=y[3:6]
    
    drdt=v

    r2=(y[0]**2+y[1]**2)**0.5
     
    dvdt_x= -G*M_disk*(y[0]/ ((r2**2+(h+(y[2]**2+b**2)**0.5)**2))**1.5)  - 2*np.pi*G*rho*pow(a,3)* \
            (y[0]/(r2*(a+r2)**2))  -  4*np.pi*G*rho1*pow(r_h,3)* ( ((r_h+r2)*np.log(1+r2/r_h) - r2)*y[0])/( r2**3 * (r_h+r2)**2)

    dvdt_y= -G*M_disk*(y[1]/ ((r2**2+(h+(y[2]**2+b**2)**0.5)**2))**1.5)  - 2*np.pi*G*rho*pow(a,3)* \
            (y[1]/(r2*(a+r2)**2))  -  4*np.pi*G*rho1*pow(r_h,3)* ( ((r_h+r2)*np.log(1+r2/r_h) - r2)*y[1])/( r2**3 * (r_h+r2)**2)

    dvdt_z= -G*M_disk*( y[2]*(h+(y[2]**2+b**2)**0.5) / ( (y[2]**2+b**2)**0.5)*(r2**2+(h+(y[2]**2+b**2)**0.5)**2)**1.5    )
    
    return [drdt[0],drdt[1],drdt[2],dvdt_x,dvdt_y,dvdt_z]


###########################Creating files

distribution=open('orbits/distribution.dat','a')
orbits1=open('orbits/orbits_general.dat','a')

ics=np.genfromtxt('ics.dat')
"""
1:j, Number
2:i, ID
3:mass1[0], initial mass 1
4:mass2[0], initial mass 2
5:mns1, neutron star mass 1
6:mns2, neutron star mass 1
7:times[0], supernova time 1
8:times1[0], supernova time 2
9:e0, initial eccentricity
10:a0,initial semimajor axis
11:e, final ecc.
12:a, final sma.
13:t, decay time
14,15,16:values[0], values[1], values[2], velocity CM first kick
17,18,19:values1[0], values1[1], values1[2], velocity CM second kick
20:z
"""

x=ics[:,0]

################## Arrays to initial-final distributions
x1=[]
y1=[]
z1=[]
R1=[]
logr1=[]
x2=[] 
y2=[]
z2=[]
R2=[] 
logr2=[]


for i in range(260,6189):
    
    print "integrating binary orbit #%d"%i
    
    x=ics[i,0]
    y=ics[i,1]
    z=ics[i,2]
    vx=ics[i,3]
    vy=ics[i,4]
    vz=ics[i,5]
    tsn1=ics[i,6]
    tsn2=ics[i,7]
    tdecay=ics[i,8]
    vk1x=ics[i,9]
    vk1y=ics[i,10]
    vk1z=ics[i,11]
    vk2x=ics[i,12]
    vk2y=ics[i,13]
    vk2z=ics[i,14]
    metal=ics[i,15]
    distance=ics[i,16]
    

    ############################Axis rotation of kick velocities

    
    #Euler angles
    inclination=np.arccos(2*random.uniform(0,1)-1)
    node=2*np.pi*random.uniform(0,1)
    periapsis=2*np.pi*random.uniform(0,1)
    
    #Sin and cos functions to (-Euler angles)
    cosi=np.cos(-inclination)
    sini=np.sin(-inclination)
    cosnode=np.cos(-node)
    sinnode=np.sin(-node)
    cosp=np.cos(-periapsis)
    sinp=np.sin(-periapsis)
    
    
    ########################Matrix rotation aplication first kick
    
    vxg1= (cosnode*cosi*cosp - sinnode*sinp)*vk1x + (cosnode*cosi*sinp + sinnode*cosp)*vk1y + (-cosnode*sini)*vk1z
    vyg1= (-sinnode*cosi*cosp - cosnode*sinp)*vk1x + (-sinnode*cosi*sinp + cosnode*cosp)*vk1y + (sinnode*sini)*vk1z
    vzg1= (sini*cosp)*vk1x + (sini*sinp)*vk1y + cosi*vk1z
    
    
    ########################Matrix rotation aplication second kick
    
    vxg2= (cosnode*cosi*cosp - sinnode*sinp)*vk2x + (cosnode*cosi*sinp + sinnode*cosp)*vk2y + (-cosnode*sini)*vk2z
    vyg2= (-sinnode*cosi*cosp - cosnode*sinp)*vk2x + (-sinnode*cosi*sinp + cosnode*cosp)*vk2y + (sinnode*sini)*vk2z
    vzg2= (sini*cosp)*vk2x + (sini*sinp)*vk2y + cosi*vk2z
    

    ###########################Differential ecuations solver
    ###        tsn1 > tsn2
    

    orbits=open('orbits/orbits_%04d.dat'%i,'a')
    times=open('orbits/times_%04d.dat'%i,'a')

    if (tsn1 > tsn2):
	## Time intervals
        t0=linspace(t,(tsn2/13)+t, 40)
        t1=linspace((tsn2/13)+t,(tsn1/13)+t,400)
        t2=linspace((tsn1/13)+t,(tdecay/13) + t,400 )
        
	# Array with initial position and velocity
        y0=([ x/3.5, y/3.5, z/3.5, vx/263, vy/263,vz/263])
        solution0=odeint(EOS,y0,t0,mxstep=400000)
        r0=solution0[:,:6]
        
        # Array with position and velocity before first kick
        y1=([r0[-1,0],r0[-1,1],r0[-1,2],r0[-1,3]+vxg1/263,r0[-1,4]+vyg1/263,r0[-1,5]+vzg1/263])
        solution1=odeint(EOS,y1,t1,mxstep=400000)
        r1=solution1[:,:6]

        # Array with position and velocity before Second kick
        y2=([r1[-1,0],r1[-1,1],r1[-1,2],r1[-1,3]+vxg2/263,r1[-1,4]+vyg2/263,r1[-1,5]+vzg2/263])
        solution2=odeint(EOS,y2,t2,mxstep=40000)
        r2=solution2[:,:6]


        ######################### Considering the maximum galactocentric distance = 50 Kpc
        r=( r2[-1,0]**2 + r2[-1,1]**2 + r2[-1,2]**2)**0.5
        print "Final galactocentric distance = %.5f"%(r*3.5)
        print "\n"
        
        ###########################Saving orbits in files
        
        x_orb=[]
        y_orb=[]
        z_orb=[]
        for j in range(len(r0[:,0])):
            x_orb.append(r0[j,0]*3.5)
            y_orb.append(r0[j,1]*3.5)
            z_orb.append(r0[j,2]*3.5)
            
            
            orbits.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            orbits1.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            
            
            if (13*t0[j]==tsn2):
                times.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))         
                
        x_orb=[]
        y_orb=[]
        z_orb=[]
        for j in range(len(r1[:,0])):
            
            x_orb.append(r1[j,0]*3.5)
            y_orb.append(r1[j,1]*3.5)
            z_orb.append(r1[j,2]*3.5)
            
            orbits.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            orbits1.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            
            
            if (13*t1[j]==tsn1):
                times.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
                
        x_orb=[]
        y_orb=[]
        z_orb=[]
        for j in range(len(r2[:,0])):
            
            x_orb.append(r2[j,0]*3.5)
            y_orb.append(r2[j,1]*3.5)
            z_orb.append(r2[j,2]*3.5)
            
            orbits.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            orbits1.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))        
            
            if (13*t2[j]==tdecay):
                times.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
                
                
        orbits.write("\n")
        orbits1.write("\n")


    else:
        ## Time intervals	
        t0=linspace(t,(tsn1/13)+t,40)
        t1=linspace((tsn1/13)+t,(tsn2/13)+t,400)
        t2=linspace((tsn2/13)+t,(tdecay/13)+t,400)
        
        # Array with initial position and velocity
        y0=([ x/3.5, y/3.5, z/3.5, vx/263, vy/263,vz/263])
        solution0=odeint(EOS,y0,t0,mxstep=400000)
        r0=solution0[:,:6]
        v0=solution0[:,3:6]
        
        # Array with position and velocity before first kick
        y1=([r0[-1,0],r0[-1,1],r0[-1,2],r0[-1,3]+vxg1/263,r0[-1,4]+vyg1/263,r0[-1,5]+vzg1/263])
        solution1=odeint(EOS,y1,t1,mxstep=400000)
        r1=solution1[:,:6]
        v1=solution1[:,3:6]
        
        # Array with position and velocity before second kick
        y2=([r1[-1,0],r1[-1,1],r1[-1,2],r1[-1,3]+vxg2/263,r1[-1,4]+vyg2/263,r1[-1,5]+vzg2/263])
        solution2=odeint(EOS,y2,t2,mxstep=40000)
        r2=solution2[:,:3]
        v2=solution2[:,3:6]
        
        #################### Considering the maximum galactocentric distance = 50 Kpc  
        r=( r2[-1,0]**2 + r2[-1,1]**2 + r2[-1,2]**2)**0.5
        print "Final galactocentric distance = %.5f"%(r*3.5)
        print "\n"
               
        ###########################Saving orbits in files
        ###     tsn1 < tsn2
        
        x_orb=[]
        y_orb=[]
        z_orb=[]
        for j in range(len(r0[:,0])):
            x_orb.append(r0[j,0]*3.5)
            y_orb.append(r0[j,1]*3.5)
            z_orb.append(r0[j,2]*3.5)
            
            orbits.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            orbits1.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            
            
            if (13*t0[j]==tsn1):
                times.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))

        x_orb=[]
        y_orb=[]
        z_orb=[]
        for j in range(len(r1[:,0])):
            
            x_orb.append(r1[j,0]*3.5)
            y_orb.append(r1[j,1]*3.5)
            z_orb.append(r1[j,2]*3.5)
            
            orbits.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            orbits1.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            
            
            if (13*t1[j]==tsn2):
                times.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
                
                
        x_orb=[]
        y_orb=[]
        z_orb=[]
        for j in range(len(r2[:,0])):
            
            x_orb.append(r2[j,0]*3.5)
            y_orb.append(r2[j,1]*3.5)
            z_orb.append(r2[j,2]*3.5)
            
            orbits.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
            orbits1.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))        
            
            if (13*t2[j]==tdecay):
                times.write("%.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j]))
                
        orbits.write("\n")
        orbits1.write("\n")


        
