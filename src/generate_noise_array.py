import sys
sys.path.append("../../src/DRD_Wind/")
import matplotlib.pyplot as plt
import numpy as np
import json
from script.main_NoiseGenerator import generateNoise

#import numba as nb
import os

def main():
    # read data
    with open('cellOutput.json') as data_file:
        cell_data = json.loads(data_file.read())
        diameter = cell_data["diameter"]
        omega = cell_data["omega"]
        numHeartBeats = cell_data["numHeartBeats"]  
        num_const_noises = cell_data["numConstNoises"]  
        radius_y = cell_data["radiusY"]
        radius_z = cell_data["radiusZ"]
        middle = cell_data["middle"]

    with open('HeartBeatSignal.json') as hbs_file:
        hbs_data = json.loads(hbs_file.read())
        hbs_data = hbs_data[0]
        velocity_values = hbs_data['v']
        area = hbs_data['a'][0]
        t_vals = hbs_data['t']

    nu = 3.3 * 10 ** -6 # m^2/s
    # Physical dx
    dx = 2 * np.sqrt(area / np.pi) / diameter #cm
    dx = dx / 100 #m
    # Physical dt
    dt = 1 / 3 * (1 / omega - 0.5) * (dx ** 2) / nu # s
    timesteps = int(numHeartBeats / dt)
    print(f"{dt/dx=}")

    num_noises = int((timesteps + 1) / num_const_noises)
    noise_field = generateNoise(int(diameter) + 1, int(diameter) + 1, num_noises)

    for idx in range(num_noises):
        normingFactor = np.amax(np.absolute(noise_field[:,:, idx,:]))
        noise_field[:,:,idx,:]   = noise_field[:,:,idx,:] / normingFactor

    for idx in range(num_noises):
        np.save(f"data/noise_field_{idx}.npy", noise_field[:,:,idx,:], allow_pickle=False)

    meta_data = {}
    meta_data["timesteps"] = timesteps
    meta_data["dt"] = dt
    meta_data["dx"] = dx
    meta_data["nu"] = nu
    meta_data["radiusY"] = radius_y
    meta_data["radiusZ"] = radius_z
    meta_data["middle"] = [middle[0], middle[1], middle[2]]

    with open("metaData.json", "w") as json_file:
        json.dump(meta_data, json_file)
    
    t_vals = [19.0] + t_vals
    t_vals = [t - 19 for t in t_vals]

    # Delete duplicate values
    t_vals = np.asarray(t_vals)
    _, unique_indices = np.unique(t_vals, return_index=True) 
    t_vals = t_vals[unique_indices]

    # Similarly, we have to add the value also to the velocity values
    # As the profile should be periodic, we simply set this value equal to the last one
    # and also delete the duplicates
    velocity_values = [velocity_values[-1]] + velocity_values
    velocity_values = np.asarray(velocity_values)
    velocity_values = velocity_values[unique_indices]

    # Convert from cm/s to m/s
    velocity_values = velocity_values / 100

    t = []
    v = []
    for idx in range(len(t_vals)):
        t.append(t_vals[idx])
        v.append(velocity_values[idx])
    print(t[0])
    print(v[0])

    with open("cleanHBData.json", "w") as json_file:
        clean = {}
        clean["t"] = t
        clean["v"] = v

        json.dump(clean, json_file)


if __name__== "__main__":
    main()