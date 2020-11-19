from amuse.units import units
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
# Order of particles: Disk => Bulge => Halo

galaxy2 = new_galactics_model(n_halo,
                                converter,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge,
                                disk_number_of_particles=n_disk)
galaxy2.move_to_center()



# move and rotate the galaxies to their appropiate initial conditions




galaxy1.rotate( 0.,np.pi/2., 0.)
galaxy2.x  += 400 | units.kpc
galaxy2.vx += -100 |units.kms
galaxy2.vy += 50 |units.kms




# Create the the gravity solvers and add the particles


converter = nbody_system.nbody_to_si(1e12|units.MSun, 100|units.kpc)
dynamics = Gadget2(converter, number_of_workers = 2)
dynamics.parameters.epsilon_squared = (100|units.parsec)**2

set1 = dynamics.particles.add_particles(galaxy1)
set2 = dynamics.particles.add_particles(galaxy2)





channel_galaxies = dynamics.particles.new_channel_to(dynamics.particles)


# The loop that evolves the gravity code

times = np.arange(0., 5001, 1) | units.Myr
threshold = 10. |units.Myr 
for time in tqdm(range(len(times))):

    dynamics.evolve_model(times[time])
    if times[time] %threshold == 0|units.Myr:
        channel_galaxies.copy()

        plt.figure(figsize = (10,8))
        plt.scatter(set1.x[:int(n_disk)].value_in(units.kpc), set1.y[:int(n_bulge)].value_in(units.kpc), s=0.3)
        plt.scatter(set1.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set1.y[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
        plt.scatter(set2.x[:int(n_disk)].value_in(units.kpc), set2.y[:int(n_bulge)].value_in(units.kpc), s=0.3)
        plt.scatter(set2.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set2.y[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
        plt.title("Disk + Bulge \n time = " +str(times[time]) )
        plt.xlabel("x [kpc]")
        plt.ylabel("y [kpc]")
        plt.xlim(-500,500)

        plt.ylim(-500,500)
        #plt.savefig("test_plots_xy/snap%04d.png"%time)
        #plt.show()
        plt.close()
        
        plt.figure(figsize = (10,8))
        plt.scatter(set1.x[:int(n_disk)].value_in(units.kpc), set1.z[:int(n_bulge)].value_in(units.kpc), s=0.3)
        plt.scatter(set1.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set1.z[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
        plt.scatter(set2.x[:int(n_disk)].value_in(units.kpc), set2.z[:int(n_bulge)].value_in(units.kpc), s=0.3)
        plt.scatter(set2.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set2.z[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
        plt.title("Disk + Bulge \n time = " +str(times[time]) )
        plt.xlabel("x [kpc]")
        plt.ylabel("z [kpc]")
        plt.xlim(-500,500)

        plt.ylim(-500,500)
        #plt.savefig("test_plots_xz/snap%04d.png"%time)
        #plt.show()
        plt.close()
        
        
        file = open("time.txt","w")
        file.write(str(times[time]))
        file.write("\n")
        #put the last data in fits files
        col1 = fits.Column(name='set_1_x', format='E', array=np.array(set1.x.value_in(units.m)))
        col2 = fits.Column(name='set_1_y', format='E', array=np.array(set1.y.value_in(units.m)))
        col3 = fits.Column(name='set_1_z', format='E', array=np.array(set1.z.value_in(units.m)))
        col4 = fits.Column(name='set_1_vx', format='E', array=np.array(set1.vx.value_in(units.kms)))
        col5 = fits.Column(name='set_1_vy', format='E', array=np.array(set1.vy.value_in(units.kms)))
        col6 = fits.Column(name='set_1_vz', format='E', array=np.array(set1.vz.value_in(units.kms)))

        
        col11 = fits.Column(name='set_2_x', format='E', array=np.array(set2.x.value_in(units.m)))
        col12 = fits.Column(name='set_2_y', format='E', array=np.array(set2.y.value_in(units.m)))
        col13 = fits.Column(name='set_2_z', format='E', array=np.array(set2.z.value_in(units.m)))
        col14 = fits.Column(name='set_2_vx', format='E', array=np.array(set2.vx.value_in(units.kms)))
        col15 = fits.Column(name='set_2_vy', format='E', array=np.array(set2.vy.value_in(units.kms)))
        col16 = fits.Column(name='set_2_vz', format='E', array=np.array(set2.vz.value_in(units.kms)))

        cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col11,col12,col13,col14,col15,col16])
        hdu = fits.BinTableHDU.from_columns(cols)
        hdu.writeto('data_large_set.fits',overwrite =True)
        file.write(str(times[time]))
        file.close()

dynamics.stop()

