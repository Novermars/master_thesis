import sys
sys.path.append("../../src/DRD_Wind/")
import matplotlib.pyplot as plt
import numpy as np
import json
from script.main_NoiseGenerator import generateNoise
from scipy.interpolate import CubicSpline
import os 

def inflow_profile(points_y, points_z, points_t, pipe_rad=1, pipe_cs_cent=[0,0,0], gamma=9):
    r_sq   = np.zeros((points_y.size,points_z.size))
    for yidx in range(len(points_y)):
        for zidx in range(len(points_z)):
            r_sq[yidx, zidx] = np.square(points_y[yidx] - pipe_cs_cent[1]) + np.square(points_z[zidx] - pipe_cs_cent[2])

    vel    = np.zeros((points_y.size,points_z.size,points_t.size,3))
    for i in range(points_t.size):
        vel[:,:,i,2] = np.maximum((gamma + 2) / gamma * (1 - (r_sq / pipe_rad ** 2) ** gamma), 0.)
    return vel

def write_json(cells_original, t_idx, turbulent_vel, indices_y, indices_z):
    data = []
    filename = f"data/inflow_profile_{t_idx}.json"
    for idx in range(len(cells_original)):
        coords = cells_original[idx]
        # Find y,z coordinate in the circle
        y = int(indices_y[coords[1]])
        z = int(indices_z[coords[2]])

        tmp = [coords[0], 
               coords[1], 
               coords[2],
               -turbulent_vel[y,z,0,2], # Inflow velocity in the +x direction
                turbulent_vel[y,z,0,0], # Inflow velocity in the +y direction
                turbulent_vel[y,z,0,1]] # Inflow velocity in the +z direction
        data.append(tmp)

    with open(filename, "w") as json_file:
        json.dump(data, json_file)
    

def main():
    # read data
    with open('cellOutput.json') as data_file:
        cell_data = json.loads(data_file.read())
        #timesteps = cell_data["timesteps"]
        diameter = cell_data["diameter"]
        omega = cell_data["omega"]
        radius = diameter / 2
        middle = cell_data["middle"]
        numHeartBeats = cell_data["numHeartBeats"]  
        num_const_noises = cell_data["numConstNoises"]  
        num_const_inflow = cell_data["numConstInflow"]
        cells_original = cell_data["cells"].copy()

    with open('HeartBeatSignal.json') as hbs_file:
        hbs_data = json.loads(hbs_file.read())
        hbs_data = hbs_data[0]
        velocity_values = hbs_data['v']
        area = hbs_data['a'][0]
        t_vals = hbs_data['t']

    print(numHeartBeats)
    # Convert to lattice units
    # Kinematic viscosity
    nu = 3.3 * 10 ** -6 # m^2/s
    print(f"area={area}")
    # Physical dx
    dx = 2 * np.sqrt(area / np.pi) / diameter #cm
    dx = dx / 100 #m
    # Physical dt
    dt = 1 / 3 * (1 / omega - 0.5) * (dx ** 2) / nu # s
    timesteps = int(numHeartBeats / dt)
    print(f"numHeartBeats={numHeartBeats}")
    print(f"dx = {dx}m")
    print(f"dt = {dt}s")
    print(timesteps)

    # Map the time values from [19.0001, 20] to [0, 1]
    # For this we have to add the value [19]
    t_vals = [19.0] + t_vals
    t_vals = [t - 19 for t in t_vals]

    # Delete duplicate values
    t_vals = np.asarray(t_vals)
    _, unique_indices = np.unique(t_vals, return_index=True) 
    t_vals = t_vals[unique_indices]
    print(len(t_vals))

    # Similarly, we have to add the value also to the velocity values
    # As the profile should be periodic, we simply set this value equal to the last one
    # and also delete the duplicates
    velocity_values = [velocity_values[-1]] + velocity_values
    velocity_values = np.asarray(velocity_values)
    velocity_values = velocity_values[unique_indices]

    # Convert from cm/s to m/s
    velocity_values = velocity_values / 100
    print(len(velocity_values))

    # Interpolate the velocity values
    interpolated_v_values = CubicSpline(t_vals, velocity_values, bc_type='periodic')
    x_vals = np.linspace(0, 1, 1000)

    # Determine the amount of noise samples to generate
    num_noises = int((timesteps + 1) / num_const_noises) + 2
    print(f"num_noises={num_noises}")
    # Generate the noise field
    noise_field = generateNoise(int(diameter) + 1, int(diameter) + 1, num_noises)

    noiseLevel    = 1./400
    normingFactor = np.amax(np.absolute(noise_field))
    noise_field   = noiseLevel * noise_field / normingFactor

    # Center the noise array cells around the middle cell given from waLBerla
    points_y = np.arange(middle[1] - radius, middle[1] + radius + 1)
    points_z = np.arange(middle[2] - radius, middle[2] + radius + 1)
    points_t = np.arange(0,1)

    # Construct hash map to map waLBerla coordinates to Python coordinates
    indices_y = {}
    indices_z = {}
    for idx, (y, z) in enumerate(zip(points_y, points_z)):
        indices_y[y] = idx
        indices_z[z] = idx

    # Construct inflow profile
    inflow_profile_data = inflow_profile(points_y, points_z, points_t, radius, middle, gamma=9)

    # In order to incorporate x and y direction noise, copy the parabolic array also to x and y
    # components + scale them if desired
    y_z_noise_factor = 1   # A value of 1 is what Ustim originally did. A value of ~10 looks really nice
                            # but is probably unphysical... no idea until we get some input from medicine
    scale_profile = np.copy(inflow_profile_data)
    scale_profile[:,:,:,0] = y_z_noise_factor*scale_profile[:,:,:,2]
    scale_profile[:,:,:,1] = y_z_noise_factor*scale_profile[:,:,:,2]

    # Put some meta data from the simulation into a .json file
    meta_data = {}
    meta_data["timesteps"] = timesteps
    meta_data["dt"] = dt
    meta_data["dx"] = dx
    meta_data["noise_factor"] = y_z_noise_factor
    meta_data["nu"] = nu

    with open("metaData.json", "w") as json_file:
        json.dump(meta_data, json_file)

    # For every timestep calculate the inflow profile and write it to a json file
    noise_idx = 0
    for t_idx in range(0, timesteps + 1 + num_const_inflow, num_const_inflow):
        velocity_value = interpolated_v_values(t_idx * dt)
        inflow_vel = velocity_value * inflow_profile_data
        cut = np.copy(noise_field[:,:,[noise_idx,noise_idx],:])
        cut = np.ascontiguousarray(cut)
        turbulent_vel = inflow_vel - velocity_value * np.multiply(cut, scale_profile)
        turbulent_vel = (dt / dx) * turbulent_vel
        write_json(cells_original, int(t_idx / num_const_inflow), turbulent_vel, indices_y, indices_z)
        
        # To save computational power, we assume that the noise is constant
        # for num_const_noises timesteps
        if t_idx % num_const_noises == 0:
            noise_idx += 1

if __name__== "__main__":
        main()