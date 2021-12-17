################################################################################
### MODULES
################################################################################

# Adapt here to your local file system's path, e.g. sys.path.append("~/.../DRD_Wind")
import sys
sys.path.append("src/DRD_Wind/")

#import os
#from numpy.core.fromnumeric import reshape
import numpy as np
import json                        # Reading time profile of heart beat in JSON format
from pyevtk.hl import gridToVTK    # Paraview output
import matplotlib.pyplot as plt



################################################################################
### MODULE FUNCTIONS
################################################################################

# Computes velocity Poiseuille profile over the square of points (points_x_i, points_y_j)_ij
# where everything outside the circle with midpoint pipe_cs_cent and radius pipe_rad is set to 0.
# The result is a Nx x Ny x Nz x 3 array where the last dimension stands for x,y and z components
# and in the third (Nz) dimension just copies of the profile are created.
# Both the x and y velocity components are zero and the copies are redundant as well. It is just made
# in this way because then one can easily add the noise (which does have x,y components and which 
# does change in z direction) and make a nice output to Paraview
def Poiseuille_profile(points_x, points_y, points_zt, pipe_rad=1, pipe_cs_cent=[0,0]):
    r_sq_x = np.square(points_x-pipe_cs_cent[0])
    r_sq_y = np.square(points_y-pipe_cs_cent[1])
    r_sq   = np.zeros((points_x.size,points_y.size))
    vel    = np.zeros((points_x.size,points_y.size,points_zt.size,3))
    for i in range(r_sq_x.size):
        r_sq[i,:] = r_sq_y + r_sq_x[i]
    for i in range(points_zt.size):
        vel[:,:,i,2] = np.maximum((pipe_rad**2 - r_sq), 0.)
    return vel


# From the cell based grid points (cell centers), create the delimiting grid by shifting
# of h/2 and appending a final point at the end (h=cell length)
def Cell2PointGrid(cellCenters):
    h = 0.5*(cellCenters[1]-cellCenters[0])
    points = cellCenters-h
    points = np.append(points,points[-1]+2.0*h)
    return points


# Output the velocity fields to Paraview (.vtr files)
def SaveToParaview(points_x, points_y, points_zt, poi_field, turb_field, tind):
    file_name = "../data/velocities"+str(tind)
    # Copy x,y and z components of velocity fields (pyevtk needs that? Or I was too
    # braindead to do it right... anyway, it works and is just for the Paraview output)
    poi_field_x_vtk = np.copy(poi_field[:,:,:,0])
    poi_field_y_vtk = np.copy(poi_field[:,:,:,1])
    poi_field_z_vtk = np.copy(poi_field[:,:,:,2])

    turb_field_x_vtk = np.copy(turb_field[:,:,:,0])
    turb_field_y_vtk = np.copy(turb_field[:,:,:,1])
    turb_field_z_vtk = np.copy(turb_field[:,:,:,2])

    cell_data = {
        "poiseuille_x": poi_field_x_vtk,
        "poiseuille_y": poi_field_y_vtk,
        "poiseuille_z": poi_field_z_vtk,

        "turbulent_x": turb_field_x_vtk,
        "turbulent_y": turb_field_y_vtk,
        "turbulent_z": turb_field_z_vtk,        
    }
    # From cell to point grid (which delimits the cells, this is definitely needed
    # by gridToVTK)
    points_x  = Cell2PointGrid(points_x)
    points_y  = Cell2PointGrid(points_y)
    points_zt = Cell2PointGrid(points_zt)
    gridToVTK(file_name, points_x, points_y, points_zt, cellData=cell_data)



################################################################################
### Main function
################################################################################
def Main():

    # Load noise field data
    # 250 x 250 x 250 x 3
    noise_field = np.load("src/DRD_Wind/data/wind.npy")

    # Set the noise level for the bloodflow simulation
    noiseLevel    = 1./4
    normingFactor = np.amax(np.absolute(noise_field))
    noise_field   = noiseLevel * noise_field/normingFactor

    # Center the noise array cells around 0 in x and y direction 
    points_x = np.arange(noise_field.shape[0])-(noise_field.shape[0]-1)/2.0
    points_y = np.arange(noise_field.shape[1])-(noise_field.shape[1]-1)/2.0
    # Scale them to the square [-1,1]^2
    points_x = points_x/points_x[-1]
    points_y = points_y/points_y[-1]
    points_zt = np.arange(0,2)


    # Load time varying velocity profile from 1D heart beat simulation
    with open('src/DRD_Wind/HeartBeatSignal.json') as hbs_file:
        data = json.loads(hbs_file.read())
        data = data[0]
        time_signal = data['v']

    # Sample time signal. Since we only have a "small" number of noise instances in the 
    # zt direction, we sample the time signal to that amount as well for a time varying 
    # noise being noise_instance_at_time_j * time_signal_at_time_j
    ts_L    = len(time_signal)
    stride  = int(np.ceil(ts_L/noise_field.shape[2]))
    ind     = slice(0,ts_L,stride)
    ts_samp = time_signal[ind]
    # Also scale the time signal to 1
    ts_samp = ts_samp/np.amax(np.absolute(ts_samp))
    ts_L    = len(ts_samp)

    # Compute parabolic profile (to prescribe Poiseuille flow and to truncate e.g. the noise
    # square field to the pipe incl. wall BCs)
    parabolic_profile = Poiseuille_profile(points_x=points_x, points_y=points_y, points_zt=points_zt)
    plt.imshow(parabolic_profile[:,:,0,2])
    plt.show()
    # In order to incorporate x and y direction noise, copy the parabolic array also to x and y
    # components + scale them if desired
    x_y_noise_factor = 0.1   # A value of 1 is what Ustim originally did. A value of ~10 looks really nice
                            # but is probably unphysical... no idea until we get some input from medicine
    scale_profile = np.copy(parabolic_profile)
    scale_profile[:,:,:,0] = x_y_noise_factor*scale_profile[:,:,:,2]
    scale_profile[:,:,:,1] = x_y_noise_factor*scale_profile[:,:,:,2]
    # Compute time instances of the velocity fields
    for tind in range(ts_L):
        # Base line flow
        poiseuille_vel = ts_samp[tind]*parabolic_profile
        # Cut out the current time steps slice from the noise field
        cut = np.copy(noise_field[:,:,[tind,tind],:])
        cut = np.ascontiguousarray(cut)
        # Turbulent flow field (Base + Noise)
        turbulent_vel = poiseuille_vel - ts_samp[tind]*np.multiply(cut,scale_profile)
        # Paraview Output    
        SaveToParaview(points_x=points_x, points_y=points_y, points_zt=points_zt, 
                       poi_field=poiseuille_vel, turb_field=turbulent_vel, tind=tind)

    # Save the turbulent velocity field to be read in by ... whatever
    np.save("../data/turbulent_flow", turbulent_vel)

################################################################################
### Entry
################################################################################
if __name__ == "__main__":
    Main()
