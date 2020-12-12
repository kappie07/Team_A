import copy
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from astropy.io import fits
import time
from os import listdir
from os.path import isfile, join

def distance(x1,x2,y1,y2,z1,z2):
	x_diff = x1-x2
	y_diff = y1-y2
	z_diff = z1-z2
	dist = (x_diff**2+ y_diff**2 + z_diff**2)**0.5
	return(dist)
    

def tail_closer_then_x (index,threshold = 0.05):
	binary_candidates = []
	for i in range(len(index)-1):
		a1 = index[i]
		a2 = index[i+1:]

		dist = distance(set_2_x[a1],set_2_x[a2],set_2_y[a1],set_2_y[a2],set_2_z[a1],set_2_z[a2])
		closest = np.where(dist <= threshold)[0]
		if len(closest) > 0:
			binary_candidates.append([i,np.where(index == a2[closest][0])[0][0]])
	return(binary_candidates)


n_halo  = 1000000
n_bulge = 500000
n_disk  = 500000
tidal_index = np.loadtxt("tidal_wave_index_4_million").astype(int)
cutoff_in_kpc = -0.4e21/3.086e19

mypath = "/data2/Roewen/AMUSE/SMA_project/Data/"
onlyfiles = sorted([f for f in listdir(mypath) if isfile(join(mypath, f))])

#for i in tqdm(range(len(onlyfiles))):
binary_candidates_all = []
for i in tqdm(range(230,301)):
    timestamp = int(onlyfiles[i][-9:-5])
    file = fits.open(mypath+onlyfiles[i])
    data1 = file[1].data
    
    set_1_x = data1.field("set_1_x")/3.086e19
    set_1_y = data1.field("set_1_y")/3.086e19
    set_1_z = data1.field("set_1_z")/3.086e19
    set_1_vx = data1.field("set_1_vx")
    set_1_vy = data1.field("set_1_vy")
    set_1_vz = data1.field("set_1_vz")
    
    set_2_x = data1.field("set_2_x")/3.086e19
    set_2_y = data1.field("set_2_y")/3.086e19
    set_2_z = data1.field("set_2_z")/3.086e19
    set_2_vx = data1.field("set_2_vx")
    set_2_vy = data1.field("set_2_vy")
    set_2_vz = data1.field("set_2_vz")
    
    binary_candidates_all.append(tail_closer_then_x(tidal_index,threshold = 0.1))
    np.savetxt("binary_candidates//bin_candidates_non_tidal_index_100pc/bin_candidates_"+str(timestamp)+".txt",binary_candidates_all[i-230],fmt = "%f")
    """
    plt.scatter(set_1_x[:n_disk],set_1_y[:n_disk],s = 0.1)
    plt.scatter(set_2_x[:n_disk],set_2_y[:n_disk],s = 0.1)
    plt.scatter(set_1_x[n_disk:n_disk+n_bulge],set_1_y[n_disk:n_disk+n_bulge],s = 0.1)
    plt.scatter(set_2_x[n_disk:n_disk+n_bulge],set_2_y[n_disk:n_disk+n_bulge],s = 0.1)
    plt.scatter(set_2_x[tidal_index],set_2_y[tidal_index],s = 0.3)
    plt.axhline(cutoff_in_kpc,linestyle = "--")    
    plt.tight_layout()
    plt.show()
	"""
