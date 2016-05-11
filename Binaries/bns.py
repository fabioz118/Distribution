#Evolve binary neutron stars using AMUSE packages and bse source codes. (Last version 11th May, 2016)
#This code need the files and fortran programs in the file "Binaries to run.

#Initial files generator with purpose to run bse.f (binary stellar evolution) code from amuse. Testing DISRUPT, both binary components Neutron Stars and decay timescale t<10000 Myr.


from amuse.lab import *
import commands
import numpy
import random
from os import system
from itertools import islice
from time import time as ptime

k=0
j=0

system("mv bns.dat bns.dat.save.$(date +%s)")
tini=ptime()
while True:
    k+=1
    
    #Constants from stellar evolution, common envelope and kick velocity
    time = 1500
    K1 = 1
    K2 = 1
    neta = 0.5
    bwind = 0
    hewind = 1
    alpha1 = 1
    lambda1 = 0.5
    ceflag = 0
    tflag = 1
    ifflag = 0
    wdflag = 1
    bhflag = 0
    nsflag = 1
    mxns = 3
    pts1 = 0.05
    pts2 = 0.01
    pts3 = 0.02
    sigma = 190.0
    beta = 0.125
    xi = 1.0
    acc2 = 1.5
    epsnov = 0.001
    eddfacc = 10.0
    gamma = -1.0

    #random seed for kick velocity
    idum = random.uniform(1,1000) 
    #if show:print "\t","Kick seed: %d",idum

    #Initial mass from Salpeter IMF using amuse packages
    mZAMS = new_salpeter_mass_distribution( 1, 8 | units.MSun, 40 | units.MSun)
    mZAMS1 = new_salpeter_mass_distribution(1, 8 |units.MSun, 40 | units.MSun)
    primary = Particles( mass = mZAMS )
    secondary = Particles( mass = mZAMS1)
    primary_mass = primary.mass.value_in( units.MSun )
    secondary_mass = secondary.mass.value_in( units.MSun)
 
    #Metalicity
    z = random.uniform(0.0004,0.03)
    
    #Eccentricity
    ecc = random.uniform( 0 , 1 )

    #Orbital period of binary system
    period = 10**random.uniform(numpy.log10(1.0),numpy.log10(10000.0)) # days

    #if show:print "\t","Kick seed: %d",idum

    files = open( "binary.in" , "w" )
    files.write( "%.5f %.5f\n %.0f\n %.5f\n %d %d\n %.5f %.5f\n %.5f %.5f %.5f\n %.5f %.5f\n %.0f %.0f %.0f %.0f %.0f %.0f %.5f\n %.0f\n %.2f %.2f %.2f\n %.5f %.5f\n %.5f %.5f %.5f %.5f %.5f "%(primary_mass, secondary_mass, time, period, K1, K2 ,z, ecc, neta, bwind,hewind,alpha1, lambda1, ceflag, tflag,ifflag,wdflag,bhflag, nsflag, mxns, idum,pts1, pts2, pts3, sigma, beta, xi, acc2, epsnov, eddfacc, gamma))
    files.close()

    """
    M = 10.017 16.029
    T = 1500
    P = 8.151
    K = 1 1
    z,e = 0.000 0.688
    wind = 0.500 0.000 1.000
    accr = 1.000 0.500
    accr = 0 1 0 1 0 1 3.000
    seed = 645
    sequence =  0.05 0.01 0.02
    sigma,beta =  190.000 0.125
    param = 1.000 1.500 0.001 10.000 -1.000     
    """

    
    #Run bse
    commands.getoutput('./bse > stars.log')

    
    #Read datafile binary.dat, which is the bse's output
     data = numpy.genfromtxt("binary.dat", dtype = 'float',invalid_raise = False)

    
    time = data[:,0]
    stellar_type1 = data[:,1]
    stellar_type2 = data[:,2]
    mass1 = data[:,3]
    mass2 = data[:,4]
    cmass1 = data[:,5]
    cmass2 = data[:,6]
    logL1 = data[:,7]
    logL2 = data[:,8]
    logr1 = data[:,9]
    logr2 = data[:,10]
    logT1 = data[:,11]
    logT2 = data[:,12]
    abin = data[:,17]
    ebin = data[:,18]

    st1 = stellar_type1
    st2 = stellar_type2
    condition = (st1==13)
    condition1 = (st2==13)
    times = time[condition]
    times1 = time[condition1]
    mns1 = cmass1[-1]
    mns2 = cmass2[-1]
    e0 = ebin[0]
    a0 = abin[0]
    e = ebin[-1]
    a = abin[-1]

    
    #id number for the binary system. 
    i = random.randint( 0,1000000 ) 
    
    #Verifying disruption via eccentricity: Values <0 indicates disruption of the system
    if e>0:
        tend=ptime()
        tela=tend-tini
        print "\t"*0,"Generating test binary %d / selected %d..."%(k,j)
        print "\t"*1,"Elapsed time: %.2f minutes, %.2f hours..."%((tend-tini)/60.0,(tend-tini)/3600.0)
        if j>0:
            irate=((tfind-tini)/3600.0)/j
        else:irate=0
        srate=(1.0*k)/(tela/3600.0)
        print "\t"*1,"Finding rate: %.2f hour/nsb..."%(irate)
        print "\t"*1,"Simulation rate: %.2f ns/hour..."%(srate)
        print "\t"*2,"Eccentricity test... Done."
        
        #Verifying stellar type for both componentes of the system

        if ( st1[-1]==13 and st2[-1]== 13):
            print "\t"*2,"Neutron star binary... Done."

            #Decay timescale by gravitational waves emission in Myr
            C = 1.3e2
            t = C * (((a**4)/(mns1*mns2*(mns1+mns2))) * pow((1-e**2),3.5)) / (1 + 3.04*e**2 + 0.4*e**4)  
            

            if t<10000:
                print "\t"*2,"Decay time test... Done."

                commands.getoutput("cp binary.in bns/binary_%06d.in"%i)
                commands.getoutput("cp binary.dat bns/binary_%06d.dat"%i)
                commands.getoutput("cp stars.log bns/binary_%06d.log"%i)

                fbin=open("binaries.dat","a")
                with open("stars.log") as myfile:

                    head = list(islice(myfile, 3))
                    values = head[1].split()
                    values1 = head[2].split()
                    fbin.write("%-6d %-08d %-15.5f %-15.5f %-15.5f %-15.5f %-15.5f %-15.5f  %-15.5f %-15.5f %-15.5f %-15.5f %-15.5f %-20s %-20s %-20s %-20s %-20s %-20s %-15.5f\n"%(j, i,mass1[0], mass2[0], mns1, mns2, times[0], times1[0], e0, a0, e, a, t, values[0], values[1], values[2], values1[0], values1[1], values1[2],z))
                j+=1
                tfind=ptime()

    
"""
Fields:    
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
