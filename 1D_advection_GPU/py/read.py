import numpy as np
import matplotlib.pyplot as plt
import json

data_dir = '../data/'
 
with open(data_dir + 'grid.json', 'r') as f:
    grid = json.load(f)
    
x = np.fromfile(data_dir + 'grid.bin', dtype=np.float32)

plt.close('all')
fig = plt.figure('advection', figsize=(5, 5))

for n in range(0,11):
    ax = fig.add_subplot(111)
    qq = np.fromfile(data_dir + 'quantity.'+str(n).zfill(8)+'.bin', dtype=np.float32)
    
    ax.plot(x, qq)
    ax.set_ylim(0,1.1)
    plt.pause(0.1)
    
    if n != 10:
        plt.clf()
    
    
    