from amuse.units import units
from amuse.units import constants
from amuse.lab import Huayno, nbody_system, new_galactics_model
import copy
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from amuse.couple import bridge
from amuse.lab import Gadget2
from amuse.community.hermite.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from astropy.io import fits
import time
from amuse.lab import Particles
from os import listdir
from os.path import isfile, join


n_halo  = 100000
n_bulge = 50000
n_disk  = 50000

def binaries(particles,hardness = 10,G = constants.G,stars_index = int(n_disk)):
    # the number over which we want to check if they're binaries
    #n = len(particles)
    n = stars_index
    #n = 10000
    
    total_particle_mass = (particles.mass).sum()
    #print(particles.mass.value_in(units.MSun),(particles.vx**2+particles.vy**2+particles.vz**2).value_in(units.kms**2))
    total_Ek = (0.5*particles.mass*(particles.vx**2+particles.vy**2+particles.vz**2)).sum()
    #print(total_Ek.value_in(units.kms**2*units.kg))
    average_Ek = total_Ek/(particles.mass).sum()
    max_mass = particles.mass.amax()
    limitE = hardness*average_Ek
    a = np.argsort(particles[int(len(particles)/2):int(len(particles)/2+n)].x.number)
    binaries = []
    #print(limitE.value_in(units.kms**2))
    #for i in tqdm(range(n-1)):
    for i in range(n):
        j = i+1
        eb_max = 0 |units.kms**2
        while j < n and (particles.x[a[j]]-particles.x[a[i]])<2*G*max_mass/limitE:
            r2=(particles.x[a[j]]-particles.x[a[i]])**2+\
               (particles.y[a[j]]-particles.y[a[i]])**2+\
               (particles.z[a[j]]-particles.z[a[i]])**2
            v2=(particles.vx[a[j]]-particles.vx[a[i]])**2+\
               (particles.vy[a[j]]-particles.vy[a[i]])**2+\
               (particles.vz[a[j]]-particles.vz[a[i]])**2
            r = r2**0.5
            #eb = G*(particles.mass[a[i]]+particles.mass[a[j]])/r-0.5*v2
            mass = 1e6 |units.MSun
            #print((G*(2*mass)/r).value_in(units.kms**2), v2.value_in(units.kms**2))
            eb = G*(2*mass)/r - 0.5*v2
            #print(eb.value_in(units.kms**2),limitE.value_in(units.kms**2))
            if r.value_in(units.kpc) < 0.1:
            	print(r.value_in(units.kpc), v2.value_in(units.kms**2))
            	print(eb.value_in(units.kms**2),limitE.value_in(units.kms**2))
            if eb > eb_max:
                eb_max = eb 
            if eb > limitE:
                print("what")
                binary = particles[a[i],a[j]].copy()
                binary.hardness = eb/average_Ek
                binaries.append(binary)
            j+=1
            #if eb_max > 0 |units.kms**2:
            #	print(eb_max.value_in(units.kms**2))
           	# 	print(particles[a[i]].position.value_in(units.kpc), particles[a[j]].position.value_in(units.kpc))
            #	print(particles[a[i]].velocity.value_in(units.kms), particles[a[j]].velocity.value_in(units.kms))
            #	print()
    return(binaries)
   
"""   
n_halo  = 10000
n_bulge = 5000
n_disk  = 5000
M_galaxy = 1e12 | units.MSun
R_galaxy = 80  | units.kpc
converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)

galaxy1 = new_galactics_model(n_halo,
                                converter,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge,
                                disk_number_of_particles=n_disk)
galaxy1.move_to_center()
start_time = time.time()
binaries(galaxy1,stars_index = 2000)
print("--- %s seconds ---" % (time.time() - start_time))
#start_time = time.time()
#galaxy1.get_binaries()   
#print("--- %s seconds ---" % (time.time() - start_time))
"""

mypath = "/home/kasper/Team_A/data/"
onlyfiles = sorted([f for f in listdir(mypath) if isfile(join(mypath, f))])

#for i in tqdm(range(len(onlyfiles))):
for i in range(1):
    file = fits.open(mypath+onlyfiles[i])
    data1 = file[1].data
    
    set_1_x = data1.field("set_1_x")
    set_1_y = data1.field("set_1_y")
    set_1_z = data1.field("set_1_z")
    set_1_vx = data1.field("set_1_vx")
    set_1_vy = data1.field("set_1_vy")
    set_1_vz = data1.field("set_1_vz")
    
    set_2_x = data1.field("set_2_x")
    set_2_y = data1.field("set_2_y")
    set_2_z = data1.field("set_2_z")
    set_2_vx = data1.field("set_2_vx")
    set_2_vy = data1.field("set_2_vy")
    set_2_vz = data1.field("set_2_vz")
	
    # galaxy 1 and galaxy 2
    set_1 = Particles(len(data1))
    set_2 = Particles(len(data1))
    
	# assign the positions and velocities to the 2 galaxies
    set_1.x = set_1_x |units.m
    set_1.y = set_1_y |units.m
    set_1.z = set_1_z |units.m
    set_1.vx = set_1_vx |units.kms
    set_1.vy = set_1_vy |units.kms
    set_1.vz = set_1_vz |units.kms

    set_2.x = set_2_x |units.m
    set_2.y = set_2_y |units.m
    set_2.z = set_2_z |units.m
    set_2.vx = set_2_vx |units.kms
    set_2.vy = set_2_vy |units.kms
    set_2.vz = set_2_vz |units.kms

	# assign the mass to their respective component

    set_1[:n_disk].mass               = 8.60749885164e+35 |units.kg
    set_1[n_disk:n_disk+n_bulge].mass = 4.73984840647e+35 |units.kg
    set_1[n_disk+n_bulge:].mass       = 1.92218326371e+37 |units.kg

    set_2[:n_disk].mass               = 8.60749885164e+35 |units.kg
    set_2[n_disk:n_disk+n_bulge].mass = 4.73984840647e+35 |units.kg
    set_2[n_disk+n_bulge:].mass       = 1.92218326371e+37 |units.kg

    set_3 = Particles()
    set_3.add_particles(set_1)
    set_3.add_particles(set_2)

    binary_list = binaries(set_3)
    if len(binary_list) > 0:
    	print(binary_list)
	
    """	
    plt.clf()
    plt.scatter(set_1_x[:n_disk],set_1_y[:n_disk],s = 0.1)
    plt.scatter(set_2_x[:n_disk],set_2_y[:n_disk],s = 0.1)
    plt.scatter(set_1_x[n_disk:n_disk+n_bulge],set_1_y[n_disk:n_disk+n_bulge],s = 0.1)
    plt.scatter(set_2_x[n_disk:n_disk+n_bulge],set_2_y[n_disk:n_disk+n_bulge],s = 0.1)
    plt.tight_layout()
    plt.show()
    #plt.draw()
    #plt.pause(0.001)
    """
#print(len(np.where(set_2_y[:n_disk] <= -0.4e21)[0]))    
