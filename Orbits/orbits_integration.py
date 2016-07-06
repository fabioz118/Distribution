########### This code integrate the trajectory of binary systems under the Milky Way's Potential.
#Modified 6th July


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
#Time units : 13 Myr

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
## Miyyamoto-Nagai parameters fit.
a1=1.22151486
b1= 0.2508

#### Definition of the  system of ordinary differential equiations. The solution are positions and velocities as function of time
def EOS(y,t):
    
    r=y[0:3]
    v=y[3:6]
    
    drdt=v

    r2=(y[0]**2+y[1]**2)**0.5
     
    dvdt_x= -G*M_disk*(y[0]/ ((r2**2+(a1+(y[2]**2+b1**2)**0.5)**2))**1.5)  - 2*np.pi*G*rho*pow(a,3)* \
            (y[0]/(r2*(a+r2)**2))  -  4*np.pi*G*rho1*pow(r_h,3)* ( ((r_h+r2)*np.log(1+r2/r_h) - r2)*y[0])/( r2**3 * (r_h+r2)**2)

    dvdt_y= -G*M_disk*(y[1]/ ((r2**2+(a1+(y[2]**2+b1**2)**0.5)**2))**1.5)  - 2*np.pi*G*rho*pow(a,3)* \
            (y[1]/(r2*(a+r2)**2))  -  4*np.pi*G*rho1*pow(r_h,3)* ( ((r_h+r2)*np.log(1+r2/r_h) - r2)*y[1])/( r2**3 * (r_h+r2)**2)

    dvdt_z= -G*M_disk*( y[2]*(a1+(y[2]**2+b1**2)**0.5) / ( (y[2]**2+b**2)**0.5)*(r2**2+(h+(y[2]**2+b**2)**0.5)**2)**1.5    )
    
    return [drdt[0],drdt[1],drdt[2],dvdt_x,dvdt_y,dvdt_z]


###########################Creating files

orbits1=open('orbits/orbits_general.dat','a')

ics=np.genfromtxt('ics_disk.dat')
"""
0:i, ID
1:mass1[0], initial mass 1
2:mass2[0], initial mass 2
3:mns1, neutron star mass 1
4:mns2, neutron star mass 1
5,6,7=x,y,z
8,9,10=vx,vy,vz
11: tsn1
12:tsn2
13: t_decay
14,15,16= kick1 vx,vy,vz
17,18,19= kick2 vx,vy,vz
20:metallicity
21: distance
22: SFR

"""

i=ics[:,0]

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
    
    iden=ics[i,0]
    mass1=ics[i,1]
    mass2=ics[i,2]
    mns1=ics[i,3]
    mns2=ics[i,4]
    x=ics[i,5]
    y=ics[i,6]
    z=ics[i,7]
    vx=ics[i,8]
    vy=ics[i,9]
    vz=ics[i,10]
    tsn1=ics[i,11]
    tsn2=ics[i,12]
    tdecay=ics[i,13]
    vk1x=ics[i,14]
    vk1y=ics[i,15]
    vk1z=ics[i,16]
    vk2x=ics[i,17]
    vk2y=ics[i,18]
    vk2z=ics[i,19]
    metal=ics[i,20]
    distance=ics[i,21]
    
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
    
    #######################Initial time for integration (bigger than the minimum time to the first merger to occur 3e7 yr)
    ## Fiducial model van der Voort et al. 2014
    t=random.uniform(2.3,769.23)  #internal time units

    ###########################Differential ecuations solver
    ###        tsn1 > tsn2
    

    orbits=open('orbits/orbits_%04d.dat'%i,'a')
    final=open('final_disk.dat','a')

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
        r=( r2[-1,0]**2 + r2[-1,1]**2 )**0.5
        print "Final galactocentric distance = %.5f"%(r*3.5)
        print "\n"
        
        ###########################Saving orbits in files. Coordinates in kpc and time in Myr.
        
        x_orb=[]
        y_orb=[]
        z_orb=[]
        vx_orb=[]
        vy_orb=[]
        vz_orb=[] 
        t_orb=[]
        for j in range(len(r0[:,0])):
            x_orb.append(r0[j,0]*3.5)
            y_orb.append(r0[j,1]*3.5)
            z_orb.append(r0[j,2]*3.5)
            vx_orb.append(r0[j,3]*263)
            vy_orb.append(r0[j,4]*263)
            vz_orb.append(r0[j,5]*263)
            t_orb.append(t0[j]*13) #Time in Myr
            
            
            orbits.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j],vx_orb[j],vy_orb[j],vz_orb[j],t_orb[j]))
            
                 
                
        x_orb=[]
        y_orb=[]
        z_orb=[]
        vx_orb=[]
        vy_orb=[]
        vz_orb=[]
        t_orb=[]
        for j in range(len(r1[:,0])):
            
            x_orb.append(r1[j,0]*3.5)
            y_orb.append(r1[j,1]*3.5)
            z_orb.append(r1[j,2]*3.5)
            vx_orb.append(r1[j,3]*263)
            vy_orb.append(r1[j,4]*263)
            vz_orb.append(r1[j,5]*263)
            t_orb.append(t1[j]*13) #Time in Myr

            orbits.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j],vx_orb[j],vy_orb[j],vz_orb[j], t_orb[j]))

                
        x_orb=[]
        y_orb=[]
        z_orb=[]
        vx_orb=[]
        vy_orb=[]
        vz_orb=[]
        t_orb=[]
        for j in range(len(r2[:,0])):
            
            x_orb.append(r2[j,0]*3.5)
            y_orb.append(r2[j,1]*3.5)
            z_orb.append(r2[j,2]*3.5)
            vx_orb.append(r2[j,3]*263)
            vy_orb.append(r2[j,4]*263)
            vz_orb.append(r2[j,5]*263)
            t_orb.append(t2[j]*13) #Time in Myr

            orbits.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j],vx_orb[j],vy_orb[j],vz_orb[j], t_orb[j]))

                
        orbits.write("\n")


    else:
        ## Time intervals	
        t0=linspace(t,(tsn1/13)+t,40)
        t1=linspace((tsn1/13)+t,(tsn2/13)+t,400)
        t2=linspace((tsn2/13)+t,(tdecay/13)+t,400)
        
        # Array with initial position and velocity
        y0=([ x/3.5, y/3.5, z/3.5, vx/263, vy/263,vz/263])
        solution0=odeint(EOS,y0,t0,mxstep=400000)
        r0=solution0[:,:6]

        
        # Array with position and velocity before first kick
        y1=([r0[-1,0],r0[-1,1],r0[-1,2],r0[-1,3]+vxg1/263,r0[-1,4]+vyg1/263,r0[-1,5]+vzg1/263])
        solution1=odeint(EOS,y1,t1,mxstep=400000)
        r1=solution1[:,:6]

        
        # Array with position and velocity before second kick
        y2=([r1[-1,0],r1[-1,1],r1[-1,2],r1[-1,3]+vxg2/263,r1[-1,4]+vyg2/263,r1[-1,5]+vzg2/263])
        solution2=odeint(EOS,y2,t2,mxstep=40000)
        r2=solution2[:,:6]

        
        #################### Considering the maximum galactocentric distance = 50 Kpc  
        r=( r2[-1,0]**2 + r2[-1,1]**2)**0.5
        print "Final galactocentric distance = %.5f"%(r*3.5)
        print "\n"
               
        ###########################Saving orbits in files
        ###     tsn1 < tsn2
             
        
        x_orb=[]
        y_orb=[]
        z_orb=[]
        vx_orb=[]
        vy_orb=[]
        vz_orb=[]
        t_orb=[]
        for j in range(len(r0[:,0])):
            x_orb.append(r0[j,0]*3.5)
            y_orb.append(r0[j,1]*3.5)
            z_orb.append(r0[j,2]*3.5)
            vx_orb.append(r0[j,3]*3.5)
            vy_orb.append(r0[j,4]*3.5)
            vz_orb.append(r0[j,5]*3.5)
            t_orb.append(t0[j]*13) #Time in Myr
            
            
            orbits.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j],vx_orb[j],vy_orb[j],vz_orb[j],t_orb[j]))
            
                 
                
        x_orb=[]
        y_orb=[]
        z_orb=[]
        vx_orb=[]
        vy_orb=[]
        vz_orb=[]
        t_orb=[]
        for j in range(len(r1[:,0])):
            
            x_orb.append(r1[j,0]*3.5)
            y_orb.append(r1[j,1]*3.5)
            z_orb.append(r1[j,2]*3.5)
            vx_orb.append(r1[j,3]*3.5)
            vy_orb.append(r1[j,4]*3.5)
            vz_orb.append(r1[j,5]*3.5)
            t_orb.append(t1[j]*13) #Time in Myr

            orbits.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j],vx_orb[j],vy_orb[j],vz_orb[j], t_orb[j]))

                
        x_orb=[]
        y_orb=[]
        z_orb=[]
        vx_orb=[]
        vy_orb=[]
        vz_orb=[]
        t_orb=[]
        for j in range(len(r2[:,0])):
            
            x_orb.append(r2[j,0]*3.5)
            y_orb.append(r2[j,1]*3.5)
            z_orb.append(r2[j,2]*3.5)
            vx_orb.append(r2[j,3]*3.5)
            vy_orb.append(r2[j,4]*3.5)
            vz_orb.append(r2[j,5]*3.5)
            t_orb.append(t2[j]*13) #Time in Myr

            orbits.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n"%(x_orb[j],y_orb[j],z_orb[j],vx_orb[j],vy_orb[j],vz_orb[j], t_orb[j]))


        orbits.write("\n")

    
    r_init = (r0[0,0]**2 + r0[0,1]**2)**0.5
    r_final = (r2[-1,0]**2 + r2[-1,1]**2)**0.5
    final.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n"%(iden, mass1, mass2, mns1, mns2, r2[-1,0], r2[-1,1], r2[-1,2], r2[-1,3], r2[-1,4], r2[-1,5], r_init, r_final))
