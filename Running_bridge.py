from amuse.units import units
from amuse.lab import nbody_system, new_galactics_model
import copy
import numpy as np
import matplotlib.pyplot as plt
#from tqdm import tqdm
from amuse.couple import bridge
from amuse.community.bhtree.interface import BHTree
from astropy.io import fits


# Milky Way := Mdisk = 4.5 Mbulge | Mhalo = 100 Mbulge
# Andromeda := Mdisk = 4-7 Mbulge | Mhalo = 87 Mbulge
# Mdisk = 5 Mbulge, Mhalo = 95 Mbulge

print("everything is fine up to now")
n_halo  = 2000
n_bulge = 1000
n_disk  = 1000
M_galaxy = 1e12 | units.MSun
R_galaxy = 80  | units.kpc
converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
print("trying to create the galaxy")

galaxy1 = new_galactics_model(n_halo,
                                converter,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge,
                                disk_number_of_particles=n_disk)
# Order of particles: Disk => Bulge => Halo

M_galaxy_2 = 1e12 | units.MSun
R_galaxy_2 = 80  | units.kpc
converter_2 = nbody_system.nbody_to_si(M_galaxy_2, R_galaxy_2)
galaxy2 = new_galactics_model(n_halo,
                                converter_2,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge,
                                disk_number_of_particles=n_disk)
                                
galaxy1.rotate(0., np.pi/2, 0.)
galaxy2.x  += 400 | units.kpc
galaxy2.y  += 100 | units.kpc
galaxy2.vx += -180 |units.kms




converter = nbody_system.nbody_to_si(1.e12|units.MSun, 100|units.kpc)
dynamics = BHTree(converter,number_of_workers=1) # ph4 Does the trick, but is kinda slow
dynamics.parameters.epsilon_squared = (100|units.parsec)**2
dynamics.parameters.timestep = 1 |units.Myr
set1 = dynamics.particles.add_particles(galaxy1)
set2 = dynamics.particles.add_particles(galaxy2)


n_halo_test  = 2000
n_disk_test  = 1000
n_bulge_test = 1000

test_particles_1 = new_galactics_model(n_halo_test,
                                converter,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge_test,
                                disk_number_of_particles=n_disk_test)

test_particles_2 = new_galactics_model(n_halo_test,
                                converter_2,
                                do_scale=True,
                                bulge_number_of_particles=n_bulge_test,
                                disk_number_of_particles=n_disk_test)
                                
                                
#set the mass of the test particles to 0
test_particles_1.mass = 0 |units.kg
test_particles_2.mass = 0 |units.kg

#we're only interested in the particles from the stars, and not the halo.
test_particles_1_stars = test_particles_1[:int(n_bulge_test+n_disk_test)]
test_particles_2_stars = test_particles_2[:int(n_bulge_test+n_disk_test)]

#rotate and translate in the same way as the galaxies
test_particles_1_stars.rotate(0., np.pi/2, 0.)
test_particles_2_stars.x  += 400 | units.kpc
test_particles_2_stars.vx += -100 |units.kms
test_particles_2_stars.vy += +10 |units.kms

star_dynamics = BHTree(converter,number_of_workers=1)
star_dynamics.parameters.timestep = 1|units.Myr
star_set_1 = star_dynamics.particles.add_particles(test_particles_1_stars)
star_set_2 = star_dynamics.particles.add_particles(test_particles_2_stars)



gravity = bridge.Bridge(use_threading=False)
gravity.add_system(star_dynamics, (dynamics,) )
gravity.add_system(dynamics)
channel = star_dynamics.particles.new_channel_to(star_dynamics.particles)

gravity.timestep = 1|units.Myr
threshold = 20. |units.Myr 

times = np.arange(0., 7000, 10) | units.Myr
for time in range(len(times)):
    gravity.evolve_model(times[time])
    if times[time] %threshold == 0|units.Myr:
        channel.copy()
        #saving the orbit of the stars with mass
        plt.figure(figsize = (10,8))
        plt.scatter(set1.x[:int(n_bulge+n_disk)].value_in(units.kpc), set1.y[:int(n_bulge+n_disk)].value_in(units.kpc), s=0.3,label = 'Big galaxy')
        plt.scatter(set2.x[:int(n_bulge+n_disk)].value_in(units.kpc), set2.y[:int(n_bulge+n_disk)].value_in(units.kpc), s=0.3,label = 'Small galaxy')
        plt.title("Disk + Bulge \n time = " +str(times[time]) )
        plt.legend()
        plt.xlabel("x [kpc]")
        plt.ylabel("y [kpc]")
        plt.xlim(-500,500)
        #plt.axis("equal")
        plt.ylim(-500,500)
        plt.savefig('./merge_plots/snap%04d.png'%time)
        #plt.show()
        plt.close()
        
        #saving the orbit of the massless tracers
        plt.figure(figsize = (10,8))
        plt.scatter(star_set_1.x[:int(n_disk_test)].value_in(units.kpc), star_set_1.y[:int(n_disk_test)].value_in(units.kpc), s=1,label = 'Big galaxy (disk)',color = 'b')
        plt.scatter(star_set_2.x[:int(n_disk_test)].value_in(units.kpc), star_set_2.y[:int(n_disk_test)].value_in(units.kpc), s=1,label = 'Small galaxy (disk)', color ='k')
        plt.scatter(star_set_1.x[int(n_disk_test):int(n_bulge_test+n_disk_test)].value_in(units.kpc), star_set_1.y[int(n_disk_test):int(n_bulge_test+n_disk_test)].value_in(units.kpc), s=1,label = 'Big galaxy (bulge)',color = 'r')
        plt.scatter(star_set_2.x[int(n_disk_test):int(n_bulge_test+n_disk_test)].value_in(units.kpc), star_set_2.y[int(n_disk_test):int(n_bulge_test+n_disk_test)].value_in(units.kpc), s=1,label = 'Small galaxy bulge',color = 'magenta')
        plt.title("Disk + Bulge for test particles \n time = " +str(times[time]) )
        plt.legend()
        plt.xlabel("x [kpc]")
        plt.ylabel("y [kpc]")
        plt.xlim(-500,500)
        #plt.axis("equal")
        plt.ylim(-500,500)
        plt.savefig('./star_plots/snapshots/snap%04d.png'%time)
        plt.close()
        '''
        file = open("time.txt","w")
        file.write(str(times[time]))
        file.write("\n")
        #put the last data in fits files
        col1 = fits.Column(name='star_1_x', format='E', array=np.array(star_set_1.x.value_in(units.m)))
        col2 = fits.Column(name='star_1_y', format='E', array=np.array(star_set_1.y.value_in(units.m)))
        col3 = fits.Column(name='star_1_z', format='E', array=np.array(star_set_1.z.value_in(units.m)))
        col4 = fits.Column(name='star_1_vx', format='E', array=np.array(star_set_1.vx.value_in(units.kms)))
        col5 = fits.Column(name='star_1_vy', format='E', array=np.array(star_set_1.vy.value_in(units.kms)))
        col6 = fits.Column(name='star_1_vz', format='E', array=np.array(star_set_1.vz.value_in(units.kms)))

        col11 = fits.Column(name='star_2_x', format='E', array=np.array(star_set_2.x.value_in(units.m)))
        col12 = fits.Column(name='star_2_y', format='E', array=np.array(star_set_2.y.value_in(units.m)))
        col13 = fits.Column(name='star_2_z', format='E', array=np.array(star_set_2.z.value_in(units.m)))
        col14 = fits.Column(name='star_2_vx', format='E', array=np.array(star_set_2.vx.value_in(units.kms)))
        col15 = fits.Column(name='star_2_vy', format='E', array=np.array(star_set_2.vy.value_in(units.kms)))
        col16 = fits.Column(name='star_2_vz', format='E', array=np.array(star_set_2.vz.value_in(units.kms)))
        
        col21 = fits.Column(name='set_1_x', format='E', array=np.array(set1.x.value_in(units.m)))
        col22 = fits.Column(name='set_1_y', format='E', array=np.array(set1.y.value_in(units.m)))
        col23 = fits.Column(name='set_1_z', format='E', array=np.array(set1.z.value_in(units.m)))
        col24 = fits.Column(name='set_1_vx', format='E', array=np.array(set1.vx.value_in(units.kms)))
        col25 = fits.Column(name='set_1_vy', format='E', array=np.array(set1.vy.value_in(units.kms)))
        col26 = fits.Column(name='set_1_vz', format='E', array=np.array(set1.vz.value_in(units.kms)))

        
        col31 = fits.Column(name='set_2_x', format='E', array=np.array(set2.x.value_in(units.m)))
        col32 = fits.Column(name='set_2_y', format='E', array=np.array(set2.y.value_in(units.m)))
        col33 = fits.Column(name='set_2_z', format='E', array=np.array(set2.z.value_in(units.m)))
        col34 = fits.Column(name='set_2_vx', format='E', array=np.array(set2.vx.value_in(units.kms)))
        col35 = fits.Column(name='set_2_vy', format='E', array=np.array(set2.vy.value_in(units.kms)))
        col36 = fits.Column(name='set_2_vz', format='E', array=np.array(set2.vz.value_in(units.kms)))

        cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col11,col12,col13,col14,col15,col16,col21,col22,col23,col24,col25,col26,col31,col32,col33,col34,col35,col36])
        hdu = fits.BinTableHDU.from_columns(cols)
        hdu.writeto('data.fits',overwrite =True)
        file.write(str(times[time]))
        file.close()
        '''
       
gravity.stop()
