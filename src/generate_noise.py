import sys
sys.path.append("../../src/DRD_Wind/")
import matplotlib.pyplot as plt
import numpy as np
import json
from script.main_NoiseGenerator import generateNoise
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
               turbulent_vel[y,z,0,2],
               turbulent_vel[y,z,0,0],
               turbulent_vel[y,z,0,1]]
        data.append(tmp)

    with open(filename, "w") as json_file:
        json.dump(data, json_file)
    

def main():
    # read data
    with open('cellOutput.json') as data_file:
        cell_data = json.loads(data_file.read())
        timesteps = cell_data["timesteps"]
        diameter = cell_data["diameter"]
        radius = diameter / 2
        middle = cell_data["middle"]
        numHeartBeats = cell_data["numHeartBeats"]    
        cells_original = cell_data["cells"].copy()

    with open('HeartBeatSignal.json') as hbs_file:
        hbs_data = json.loads(hbs_file.read())
        hbs_data = hbs_data[0]
        time_signal = hbs_data['v']

    # Amount we have to skip in the time signal
    stride = len(time_signal) / (timesteps * numHeartBeats)

    noise_field = generateNoise(int(diameter) + 1, int(diameter) + 1, int(timesteps) + 1)

    noiseLevel    = 1./4
    normingFactor = np.amax(np.absolute(noise_field))
    noise_field   = noiseLevel * noise_field/normingFactor

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

    inflow_profile_data = inflow_profile(points_y, points_z, points_t, radius, middle, 9)

    # In order to incorporate x and y direction noise, copy the parabolic array also to x and y
    # components + scale them if desired
    y_z_noise_factor = 10   # A value of 1 is what Ustim originally did. A value of ~10 looks really nice
                            # but is probably unphysical... no idea until we get some input from medicine
    scale_profile = np.copy(inflow_profile_data)
    scale_profile[:,:,:,0] = y_z_noise_factor*scale_profile[:,:,:,2]
    scale_profile[:,:,:,1] = z_z_noise_factor*scale_profile[:,:,:,2]

    # Slice the array depending on the number of heart beat cycles and number of timesteps
    ind = slice(0, len(time_signal), int(stride))
    ts_samples = time_signal[ind]

    # For every timestep calculate the inflow profile and write it to a json file
    for t_idx in range(len(ts_samples)):
        inflow_vel = ts_samples[t_idx] * inflow_profile_data
        cut = np.copy(noise_field[:,:,[t_idx,t_idx],:])
        cut = np.ascontiguousarray(cut)
        turbulent_vel = inflow_vel - ts_samples[t_idx] * np.multiply(cut, scale_profile)
        write_json(cells_original, t_idx, turbulent_vel, indices_y, indices_z)

if __name__== "__main__":
        main()