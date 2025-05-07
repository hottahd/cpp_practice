import numpy as np
import matplotlib.pyplot as plt
import json

data_dir = '../data/'
 
with open(data_dir + "time_output.txt", "r") as f:
    n_output = int(f.readline().strip())
 
with open(data_dir + 'config.json', 'r') as f:
    config = json.load(f)
        
x = np.fromfile(data_dir + 'grid.bin', dtype=np.float64)
i_size = config['grid']['i_size']
j_size = config['grid']['j_size']
k_size = config['grid']['k_size']

iskip = 1 if i_size > 1 else 0
jskip = 1 if j_size > 1 else 0
kskip = 1 if k_size > 1 else 0

i_margin = config['grid']['margin']*iskip
j_margin = config['grid']['margin']*jskip
k_margin = config['grid']['margin']*kskip

i_total = i_size + 2*i_margin
j_total = j_size + 2*j_margin
k_total = k_size + 2*k_margin

plt.close('all')
fig = plt.figure('MHD', figsize=(15, 5))

n_start = 0
n_end = n_output - 1
for n in range(n_start,n_end+1):
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    qq = np.fromfile(data_dir + 'mhd.'+str(n).zfill(8)+'.bin', dtype=np.float64).reshape((9, i_total, j_total, k_total))
    ro = qq[0, :, :, :]
    vx = qq[1, :, :, :]
    vy = qq[2, :, :, :]
    vz = qq[3, :, :, :]
    bx = qq[4, :, :, :]
    by = qq[5, :, :, :]
    bz = qq[6, :, :, :]
    ei = qq[7, :, :, :]
    
    ax1.plot(ro[:, 0, 0])
    ax2.plot(vx[:, 0, 0])
    ax3.plot(ei[:, 0, 0])
    plt.pause(0.1)
    if(n != n_end):
        plt.clf()
    
    
    