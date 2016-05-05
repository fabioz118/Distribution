This text is with the purpose to explain, describe and keep record of the way that our codes and routines were executed, the requirements and the compilation characteristics (Fabio A. Cardona & Jorge I. Zuluaga 2016)

Here the steps:

1.Neutron Star Binaries generation:

To generate and evolve the binary systems, the AMUSE packages and codes are used. For this, installation of the AMUSE is needed.

Follow the installation steps described in the link:
http://amusecode.org/doc/install/howto-install-AMUSE.html

Once installed, the code bns.py can run without any problem. Also, bns.py need the source codes in the folder, codes like kick.f, bse.f, comenv.f. This codes are written in fortran, and need to be compilated using the MAKEFILE in the same folder.

bse-manual.txt explains the values needed in the initial conditions for the running. THis conditions are generated in bns.p

This codes produce a file called "binaries.dat" in which are the fundamental parameters of each binary system to continue with the procedures. Parameters, which are:

1: Number
2: ID
3: Initial mass 1
4: Initial mass 2
5: Neutron star mass 1
6: Neutron star mass 1
7: Supernova time 1
8: Supernova time 2
9: e0, Initial eccentricity
10: a0,Initial semimajor axis
11: e, Final eccentricity
12: a, Final semimajor axis
13: t, decay timescale
14,15,16: vx1,vy1, vz1 , Center of Mass velocity at first kick
17,18,19: vx2,vy2, vz2 , Center of Mass velocity at second kick
20: z Metallicity

1.1. Plots of the selected neutron star binaries

Code plots.py plot some parameters of the resultant systems. 



2. Initial conditions and Rejection


The file "ics_disk.dat" have the positions and velocities for a milky way-like exponential disk generated with GADGET 2.


