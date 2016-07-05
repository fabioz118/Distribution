####This code makes a rejection using initial positions and velocities generated with GADGET 2. Also, asign  metallicity and Star Formation Rate (SFR) according with the Milky Way radial gradients.

from amuse.lab import *
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os,commands

commands.getoutput("rm initial_conditions.txt")

np.random.seed=1


#Data file with bns
bns = np.genfromtxt('binarias.dat')
"""
0,1:
2,3 : mass1[0], mass2[0]
4,5 :  mns1, mns2
6,7 : times[0], times1[0]
8,9 :  e0, a0
10,11 : e, a
12 : t_decay 
13,14,15 : values[0](vx_kick1), values[1](vy_kick1), values[2] (vz_kick1), 
16,17,18 : values1[0],values1[1], values1[2]
19 :  Metallicity
"""
mass1 = bns[:,2]
mass2 = bns[:,3]
mns1 =  bns[:,4]
mns2 =  bns[:,5]
t1 =  bns[:,6]
t2 =  bns[:,7]
t_decay =  bns[:,12]
kickx1 =  bns[:,13]
kicky1 =  bns[:,14]
kickz1 =  bns[:,15]
kickx2 =  bns[:,16]
kicky2 =  bns[:,17]
kickz2 =  bns[:,18]
metallicity = bns[:,19]


#Data file with initial postions(kpc) and velocities(km/s)
disk = np.genfromtxt('ics_disk.dat')
disk_id = disk[:,0]
disk_x = disk[:,1]
disk_y = disk[:,2]
disk_z = disk[:,3] 
disk_vx = disk[:,4]
disk_vy = disk[:,5]
disk_vz = disk[:,6]

#Data file of SFR radial gradient
SFR = np.genfromtxt('SFR.dat')
SFR_R = SFR[:,0]
SFR_rate = SFR[:,1]

#Metallicity gradient
metals = np.genfromtxt('metallicity_gradient.dat')
R = metals[:,0]
Z = metals[:,1]

i = 0
j = 0

x=[]
y=[]
z=[]
vx=[]
vy=[]
vz=[]
iden=[]
dist=[]
form_rate=[]

## Metallicity dispersion
disp = 0.005

#####################################################################
############ Star Formation Rate (SFR) assignment ###################
#####################################################################

qload=0
if not qload:
    
    ilim=6000
    icount=0
    for i in range(len(disk_id)):
        
        #conversion to cilindric coordinates from cartesians
        disk_R = ( disk_x[i] ** 2 + disk_y[i] ** 2) ** 0.5
        
        
        for j in range(len(SFR_R)):
            
            if ( disk_R > SFR_R[j-1] and disk_R < SFR_R[j] ) :
                
                rate = random.uniform(0,6)
                
                
                if (rate < SFR_rate[j] ) :
                    
                    x.append(disk_x[i])
                    y.append(disk_y[i])
                    z.append(disk_z[i])
                    vx.append(disk_vx[i])
                    vy.append(disk_vy[i])
                    vz.append(disk_vz[i])
                    iden.append(disk_id[i])
                    dist.append(disk_R)
                    form_rate.append(SFR_rate[j])
                    icount+=1
                    
                    if icount>ilim:break
                    
            if icount>ilim:break
                    
    np.savetxt("temp.dat",np.column_stack((x,y,z,vx,vy,vz,iden,dist)))

                    
else:
    data=np.loadtxt("temp.dat")
    i=0
    x=data[:,i];i+=1 #(disk_x[i])
    y=data[:,i];i+=1 #(disk_y[i])
    z=data[:,i];i+=1 #(disk_z[i])
    vx=data[:,i];i+=1 #(disk_vx[i])
    vy=data[:,i];i+=1 #(disk_vy[i])
    vz=data[:,i];i+=1 #(disk_vz[i])
    iden=data[:,i];i+=1 #(disk_id[i])
    dist=data[:,i];i+=1 #(disk_R)


#####################################################
############## Metallicity assignment################
#####################################################


distancia=[]
qload=0
disp = 0.0
if not qload:
    metal = []
    m1 = []
    m2 = []
    distancia = []
    for i in range(len(x)):

        #IMF (Salpeter 1995)
        mZAMS = new_salpeter_mass_distribution( 1, 8 | units.MSun, 40 | units.MSun)
        mZAMS1 = new_salpeter_mass_distribution(1, 8 |units.MSun, 40 | units.MSun)
        primary = Particles( mass = mZAMS )
        secondary = Particles( mass = mZAMS1)
        primary_mass = primary.mass.value_in( units.MSun )
        secondary_mass = secondary.mass.value_in( units.MSun)

        for j in range(len(R)):

            if (dist[i] > R[j-1] and dist[i] < R[j]):

                #from iron abundances to metallicity
                metals = 0.771 / ( 3 + (1/(0.0263*pow(10,0.95*Z[j]))))
                metals = random.uniform(metals - disp,metals + disp)
                metal.append(metals)
                m1.append(primary_mass)
                m2.append(secondary_mass)
                distancia.append(R[j])
                
    np.savetxt("temp2.dat",np.column_stack((metal,m1,m2,distancia)))
    

else:
    data=np.loadtxt("temp2.dat")
    metal=data[:,0]
    m1=data[:,1]
    m2=data[:,2]
    dist=data[:,3]
    

plt.figure(1)
plt.plot(distancia,metal,'b.')
plt.savefig("metal_grad.png")


ics=open("ics.dat","a")
isel=[]

############################################################################################################
#Calculate the parametric distance between bns generated wihth bse and Gadget particles in a "m1,m2,z space
############################################################################################################

for i in range(len(mass1)):
    
    print "Comparison: ",metallicity[i],mass1[i],mass2[i]
    
    param=[]
    
    for j in range(len(metal)):

      
        # solar metallicity = 0.017, 1600 = 40 solar masses squared (maximum mass value)
        d =( 10000*( np.log10(metallicity[i]/metal[j])**2 / (np.log10(0.017/0.0001)**2)) +((mass1[i] - m1[j])**2 / 1600 ) + \
((mass2[i] - m2[j])**2 / 1600 ) )**0.5        
        param.append(d)

    
#    parametric= np.array(param)
    isort = np.argsort(param)   
    print isort

    k=0
    iclose=isort[k]
    for iclose in isel:
        print "otra"
        k+=1
        iclose=isort[k]
        isel+=[iclose]
        print "\tClosest: ",metal[iclose],m1[iclose],m2[iclose],param[iclose]


    print "sorted:", metal[isort[0]],m1[isort[0]],m2[isort[0]]
    print isort[0]
    ics.write("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f"\
%(mass1[i],mass2[i],mns1[i],mns2[i],x[isort[0]], y[isort[0]],z[isort[0]],vx[isort[0]],vy[isort[0]],vz[isort[0]],t1[i],t2[i],t_decay[i],\
 kickx1[i],kicky1[i],kickz1[i],kickx2[i],kicky2[i],kickz2[i],metal[isort[0]],dist[isort[0]],form_rate[isort[0]]))

"""
1,2,3=x,y,z
4,5,6=vx,vy,vz
7: tsn1
8:tsn2
9: t_decay
10,11,12= kick1 vx,vy,vz
13,14,15= kick2 vx,vy,vz
16:metallicity
17: distance
18: SFR
"""
    
