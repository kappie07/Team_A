import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm
from astropy.io import fits
from os import listdir
from os.path import isfile, join

n_halo = 100000
n_disk = 50000
n_bulge = 50000


mypath = "/home/kasper/Team_A/data/"
onlyfiles = sorted([f for f in listdir(mypath) if isfile(join(mypath, f))])

#for i in tqdm(range(len(onlyfiles))):
for i in range(1):
    file = fits.open(mypath+onlyfiles[i])
    data1 = file[1].data
    
    set_2_x = np.array(data1.field("set_2_x")[:n_disk],dtype = "float64")/3.086e16
    set_2_y = np.array(data1.field("set_2_y")[:n_disk],dtype = "float64")/3.086e16
    set_2_z = np.array(data1.field("set_2_z")[:n_disk],dtype = "float64")/3.086e16
    set_2_vx = np.array(data1.field("set_2_vx")[:n_disk],dtype = "float64")
    set_2_vy = np.array(data1.field("set_2_vy")[:n_disk],dtype = "float64")
    set_2_vz = np.array(data1.field("set_2_vz")[:n_disk],dtype = "float64")
    
   
closest_distances = np.zeros(n_disk)
for i in tqdm(range(len(set_2_x))):
    a1 = i
    a2 = np.delete(np.arange(n_disk),i)
    absolute_x = (set_2_x[a1] - set_2_x[a2])**2
    absolute_y = (set_2_y[a1] - set_2_y[a2])**2
    absolute_z = (set_2_z[a1] - set_2_z[a2])**2
    dist = (absolute_x + absolute_y + absolute_z)**0.5
    closest = min(dist)
    closest_distances[i] = closest

    
    
    
bins = np.linspace(0,300,300)
np.savetxt("closest_initial_distance_400000_particles.txt",closest_distances, fmt = "%f")
plt.hist(closest_distances,bins = bins)
plt.xlabel("closest particle [pc]")
plt.ylabel("counts")
plt.savefig("report_plots/initial_distance_distribution_400000_particles")
plt.show()
	   
