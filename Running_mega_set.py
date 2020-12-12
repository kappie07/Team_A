from amuse.units import units
from amuse.lab import Gadget2
import copy
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from astropy.io import fits
import time
from amuse.lab import Particles
from amuse.lab import nbody_system

n_halo  = 1000000
n_bulge = 500000
n_disk  = 500000

start_time = time.time()
def run(saved_time):
	mypath = "/home/kasper/Team_A/Tidal_Wave/Large/Data/"
	file = "Gadget2_save_data_"+str(np.round(saved_time.astype(int)))+".fits"
	file = fits.open(mypath+file)
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

	#print("data assigned")
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
		
	set_1[:n_disk].mass               = 8.60749885164e+34 |units.kg
	set_1[n_disk:n_disk+n_bulge].mass = 4.73984840647e+34 |units.kg
	set_1[n_disk+n_bulge:].mass       = 1.92218326371e+36 |units.kg

	set_2[:n_disk].mass               = 8.60749885164e+34 |units.kg
	set_2[n_disk:n_disk+n_bulge].mass = 4.73984840647e+34 |units.kg
	set_2[n_disk+n_bulge:].mass       = 1.92218326371e+36 |units.kg


	del(set_1_x)
	del(set_1_y)
	del(set_1_z)
	del(set_1_vx)
	del(set_1_vy)
	del(set_1_vz)
	del(set_2_x)
	del(set_2_y)
	del(set_2_z)
	del(set_2_vx)
	del(set_2_vy)
	del(set_2_vz)
	#plt.scatter(set_1_x[:n_disk],set_1_y[:n_disk],s = 0.1)
	#plt.scatter(set_2_x[:n_disk],set_2_y[:n_disk],s = 0.1)
	#plt.scatter(set_1_x[n_disk:n_disk+n_bulge],set_1_y[n_disk:n_disk+n_bulge],s = 0.1)
	#plt.scatter(set_2_x[n_disk:n_disk+n_bulge],set_2_y[n_disk:n_disk+n_bulge],s = 0.1)
	#plt.show()

	converter = nbody_system.nbody_to_si(1.e12|units.MSun, 100|units.kpc)
	dynamics = Gadget2(converter, number_of_workers=2)

	dynamics.parameters.epsilon_squared = (100|units.parsec)**2
	set1 = dynamics.particles.add_particles(set_1)
	set2 = dynamics.particles.add_particles(set_2)
	del(set_1)
	del(set_2)
	channel = dynamics.particles.new_channel_to(dynamics.particles)

	times = np.arange(1, 101, 1) | units.Myr
	threshold = 10 |units.Myr 

	for time in tqdm(range(0, len(times))):
		dynamics.evolve_model(times[time])
		if times[time] %threshold == 0|units.Myr:
			channel.copy()
			current_time = saved_time + times[time].value_in(units.Myr)
			np.savetxt("current_time.txt",[current_time],fmt = "%f")
			plt.figure(figsize = (10,8))
			plt.scatter(set1.x[:int(n_disk)].value_in(units.kpc), set1.y[:int(n_bulge)].value_in(units.kpc), s=0.3)
			plt.scatter(set1.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set1.y[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
			plt.scatter(set2.x[:int(n_disk)].value_in(units.kpc), set2.y[:int(n_bulge)].value_in(units.kpc), s=0.3)
			plt.scatter(set2.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set2.y[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
			plt.title("Disk + Bulge \n time = " +str(current_time) )
			#plt.legend()
			plt.xlabel("x [kpc]")
			plt.ylabel("y [kpc]")
			plt.xlim(-500,500)
			#plt.axis("equal")
			plt.ylim(-500,500)
			plt.savefig("Tidal_Wave/Large/xy/snap%04d.png"%current_time)
			#plt.show()
			plt.close()

			plt.figure(figsize = (10,8))
			plt.scatter(set1.x[:int(n_disk)].value_in(units.kpc), set1.z[:int(n_bulge)].value_in(units.kpc), s=0.3)
			plt.scatter(set1.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set1.z[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
			plt.scatter(set2.x[:int(n_disk)].value_in(units.kpc), set2.z[:int(n_bulge)].value_in(units.kpc), s=0.3)
			plt.scatter(set2.x[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set2.z[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
			plt.title("Disk + Bulge \n time = " +str(current_time) )
			#plt.legend()
			plt.xlabel("x [kpc]")
			plt.ylabel("z [kpc]")
			plt.xlim(-500,500)
			#plt.axis("equal")
			plt.ylim(-500,500)
			plt.savefig("Tidal_Wave/Large/xz/snap%04d.png"%current_time)
			#plt.show()
			plt.close()

			plt.figure(figsize = (10,8))
			plt.scatter(set1.y[:int(n_disk)].value_in(units.kpc), set1.z[:int(n_bulge)].value_in(units.kpc), s=0.3)
			plt.scatter(set1.y[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set1.z[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
			plt.scatter(set2.y[:int(n_disk)].value_in(units.kpc), set2.z[:int(n_bulge)].value_in(units.kpc), s=0.3)
			plt.scatter(set2.y[int(n_disk):int(n_disk+n_bulge)].value_in(units.kpc), set2.z[int(n_disk):int(n_bulge+n_disk)].value_in(units.kpc), s=0.3)
			plt.title("Disk + Bulge \n time = " +str(current_time) )
			#plt.legend()
			plt.xlabel("y [kpc]")
			plt.ylabel("z [kpc]")
			plt.xlim(-500,500)
			#plt.axis("equal")
			plt.ylim(-500,500)
			plt.savefig("Tidal_Wave/Large/yz/snap%04d.png"%current_time)
			#plt.show()
			plt.close()

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

			cols = fits.ColDefs([col21,col22,col23,col24,col25,col26,col31,col32,col33,col34,col35,col36])
			hdu = fits.BinTableHDU.from_columns(cols)
			name = "Tidal_Wave/Large/Data/Gadget2_save_data_%04d.fits"%current_time
			hdu.writeto(name,overwrite =True)
			del(cols)
			#del(set1)
			#del(set2)
			del(hdu)
			del(col21,col22,col23,col24,col25,col26,col31,col32,col33,col34,col35,col36)
	dynamics.stop()
	



saved_time = np.loadtxt("current_time.txt", dtype = "float")
print(saved_time)
run(saved_time)
