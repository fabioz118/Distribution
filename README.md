This text is with the purpose to explain, describe and keep record of the way that our codes and routines were executed, the requirements and the compilation characteristics (Fabio A. Cardona & Jorge I. Zuluaga 2016)

Here the steps:

###Neutron Star Binaries generation (Binaries Folder):
- Input: All the files ".f" and ".o" in the file.


To generate and evolve the binary systems, the AMUSE packages and codes are used. For this, installation of the AMUSE (http://amusecode.org/) is needed.

Follow the installation steps described in the link:
http://amusecode.org/doc/install/howto-install-AMUSE.html

Once installed, the code bns.py can run without any problem. Also, bns.py needs the source codes in the folder "Binaries", this codes are written in fortran, and need to be compilated using the MAKEFILE in the same folder.

bse-manual.txt explains the values needed in the initial conditions for the running. THis conditions are generated in bns.p

- Output:
This codes produce a file called "binaries.dat" in which are the fundamental parameters of each binary system to continue with the procedures.

..- 1: Number
- 2: ID
- 3: Initial mass 1
- 4: Initial mass 2
- 5: Neutron star mass 1
- 6: Neutron star mass 1
- 7: Supernova time 1
- 8: Supernova time 2
- 9: e0, Initial eccentricity
- 10: a0,Initial semimajor axis
- 11: e, Final eccentricity
- 12: a, Final semimajor axis
- 13: t, decay timescale
- 14,15,16: vx1,vy1, vz1 , Center of Mass velocity at first kick
- 17,18,19: vx2,vy2, vz2 , Center of Mass velocity at second kick
- 20: z Metallicity


###Plots of the selected neutron star binaries:
Input: "binaries.dat"
Output: Scatter plots ".png"

Code plots.py do some scatter plots for parameters of the resultant systems. 


###Initial conditions and Rejection (Rejection Folder):
Input: "binaries.dat", "metallicity_gradient.dat", "SFR.dat"
Output: "ics_disk.dat"

Now, code rejection.py run a rejection of initial conditions to include the  metallicity gradient and star formation rate (SFR). The output is a file with positions and velocities called "ics.dat":

*1,2,3=x,y,z
*4,5,6=vx,vy,vz
*7: supernova time 1
*8: supernova time 2
*9: decay timescale
*10,11,12= kick1 vx,vy,vz
*13,14,15= kick2 vx,vy,vz
*16: metallicity
*17: Galactocentric distance
*18: SFR


###Orbits integration (Orbits folder):
Input:"ics_disk.dat"
Output:Folder "Orbits" containing "orbits_*dat", "times_*.dat", "orbits_general.dat"

Code "orbit_integration.py" solve the equation of motion in different time intervals. The motion is under potential for a Milky Way galaxy.
The run output in directory called "orbits" the positions and velocities in each step time for each binary system. 



