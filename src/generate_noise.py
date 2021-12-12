import numpy as np
import json
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))
#print("test")
print(dir_path)

# read data
with open('cellOutput.json') as data_file:
    data = json.loads(data_file.read())
timesteps = data["timesteps"]
cells_original = data["cells"].copy()

for time in range(timesteps + 1):
    filename = f"noise_data/noise_data_{time}.json"
    cells = []
    for idx in range(len(cells_original)):
        tmp = []
        x = np.random.uniform(-0.01, 0.01)
        y = np.random.uniform(-0.01, 0.01)
        z = np.random.uniform(-0.01, 0.01)
        tmp.append(cells_original[idx][0])
        tmp.append(cells_original[idx][1])
        tmp.append(cells_original[idx][2])
        tmp.append(x)
        tmp.append(y)
        tmp.append(z)
        cells.append(tmp)

    with open(filename, "w") as json_file:
        json.dump(cells, json_file)