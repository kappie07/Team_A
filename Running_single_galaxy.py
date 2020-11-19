from amuse.units import units
from amuse.lab import Huayno, nbody_system, new_galactics_model
import copy
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from amuse.lab import Gadget2
from astropy.io import fits
import time


# Initialize the galaxies

n_halo  = 20000
n_bulge = 10000
n_disk  = 10000
M_galaxy = 1e12 | units.MSun
R_galaxy = 80  | units.kpc
converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)

galaxy1 = new_galactics_model(n_halo,
                                converter,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge,
                                disk_number_of_particles=n_disk)
galaxy1.move_to_center()


# Create the the gravity solvers and add the particles


converter = nbody_system.nbody_to_si(1e12|units.MSun, 100|units.kpc)
dynamics = Gadget2(converter, number_of_workers = 2)
dynamics.parameters.epsilon_squared = (100|units.parsec)**2
set1 = dynamics.particles.add_particles(galaxy1)

channel_galaxies = dynamics.particles.new_channel_to(dynamics.particles)


# The loop that evolves the gravity code


times = np.arange(0., 1001, 1) | units.Myr
threshold = 10. |units.Myr 
for time in tqdm(range(len(times))):

    dynamics.evolve_model(times[time])
    if times[time] %threshold == 0|units.Myr:
        channel_galaxies.copy()
        
        plt.figure(figsize = (10,8))
        plt.scatter(set1.x[:int(n_disk)].value_in(units.kpc), set1.y[:int(n_bulge)].value_in(units.kpc), s=0.3)
        plt.scatter(set1.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set1.y[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
        plt.title("Disk + Bulge \n time = " +str(times[time]) )
        plt.xlabel("x [kpc]")
        plt.ylabel("y [kpc]")
        plt.xlim(-25,25)
        plt.ylim(-25,25)
        #plt.savefig("test_plots_xy/snap%04d.png"%time)
        plt.close()
        
        plt.figure(figsize = (10,8))
        plt.scatter(set1.x[:int(n_disk)].value_in(units.kpc), set1.z[:int(n_bulge)].value_in(units.kpc), s=0.3)
        plt.scatter(set1.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set1.z[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
        plt.title("Disk + Bulge \n time = " +str(times[time]) )
        plt.xlabel("x [kpc]")
        plt.ylabel("z [kpc]")
        plt.xlim(-25,25)
        plt.ylim(-25,25)
        #plt.savefig("test_plots_xz/snap%04d.png"%time)

        plt.close()
        
dynamics.stop()
