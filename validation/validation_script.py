import yaml
import subprocess
import os
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


# validate solver against results from ghia et. al. 
# square domain (L = H = 1) and lid driven condition with U = 1

Reynold_numbers = [10000]

# go to project home directory
os.chdir('..')

for Re in Reynold_numbers:

    # open simulation data file
    with open('simdata.yaml') as f:
	dat = yaml.load(f)

    grid_size = dat['grid']['N_x']
    dx = dat['grid']['L_x']/float(dat['grid']['N_x'])
    dat['params']['Re'] = Re
   
    # CFL based computation of dt
    dat['params']['dt']  = 0.005*(dx**2)*Re
    dat['params']['nsteps'] = int(100/dat['params']['dt'])

    # update simulation data file
    with open('simdata.yaml', 'w') as f:
         yaml.dump(dat, f)    

    print "Running simulation for Reynolds number: ", "%4d"%Re 
    # run simulation:
    sim = subprocess.Popen('bin/test', stdout=None)
    sim.communicate()

    # read the results:
    u = np.loadtxt('results/u.txt')
    v = np.loadtxt('results/v.txt')
    u = u.reshape([grid_size, grid_size])    
    v = v.reshape([grid_size, grid_size])
    
    [x, y] = np.meshgrid(np.linspace(0, 1, grid_size), np.linspace(0, 1, grid_size))
    plt.streamplot(x, y, u, v)
    plt.savefig('streamplot%d.png'%Re)
    plt.close()

    # extract the centre line velocities:
    u_centre = u[:,grid_size/2]      # u velocities about central vertical line
    v_centre = v[grid_size/2,:]      # v velocities about central horizontal line

    # make plots
    u_pts = np.loadtxt('validation/ghia/u/points.txt')
    u_ghia = np.loadtxt('validation/ghia/u/%d.txt'%Re)
    
    plt.plot(u_pts, u_ghia, 'o--k', label='Ghia et. al')
    plt.plot(np.linspace(0, 1, grid_size), u_centre, 'k', label='GPU solver')
    plt.legend(loc='best')
    plt.savefig('validation/%d.png'%Re)
    plt.close()
