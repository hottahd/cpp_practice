import numpy as np
import matplotlib.pyplot as plt
import json

data_dir = '../data/'
 
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

fig = plt.figure('MHD', figsize=(5, 5))

for n in range(0,11):
    ax = fig.add_subplot(111)
    qq = np.fromfile(data_dir + 'mhd.'+str(n).zfill(8)+'.bin', dtype=np.float64).reshape((9, i_total, j_total, k_total))
    ro = qq[0, :, :, :]
    vx = qq[1, :, :, :]
    vy = qq[2, :, :, :]
    vz = qq[3, :, :, :]
    bx = qq[4, :, :, :]
    by = qq[5, :, :, :]
    bz = qq[6, :, :, :]
    ei = qq[7, :, :, :]
    
    ax.plot(ro[:, 0, 0])
    plt.pause(0.1)
    if(n != 10):
        plt.clf()
    
    
    