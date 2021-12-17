
import torch
import matplotlib.pyplot as plt
from pylab import *
import pickle


from source.Calibration import CalibrationProblem
from source.DataGenerator import OnePointSpectraDataGenerator

####################################
### Configuration
####################################

case = 'tauNet_Kaimal_pen_reg' ### 'tauNet_Kaimal_noreg', 'tauNet_Kaimal_noisy', 'tauNet_Kaimal_pen_reg', 'Mann_Kaimal', 'tauNet_Kaimal', 'tauNet_Kaimal_cross', 'tauNet_SimuYeo'
input_folder  = '/path/to/code/data/'
output_folder = '/path/to/figures/' + case + '/'


####################################
### Import results
####################################
with open(input_folder+case+'.pkl', 'rb') as file:
    config, parameters, Data, _, loss_history = pickle.load(file)

pb = CalibrationProblem(**config)
pb.parameters = parameters

print(pb.OPS.update_scales())


# ####################################
# ### Data
# ####################################
# k1_data_pts = np.logspace(-1, 2, 20)
# DataPoints  = [ (k1, 1) for k1 in k1_data_pts ]
# Data = OnePointSpectraDataGenerator(DataPoints=DataPoints, **config).Data


# ####################################
# ### Plot
# ####################################
# pb.plot(Data=Data)


####################################
### Figure data
####################################
x = np.logspace(-1, 2, 20) ### NOTE: Experiment 1: np.logspace(-1, 2, 20), Experiment 2: np.logspace(-2, 2, 40)
y = pb.eval(x)
plt.plot(x,y[0])
plt.plot(x,y[1])
plt.plot(x,y[2])
plt.plot(x,-y[3])
# plt.plot(x,y[4])
# plt.plot(x,y[5])
plt.plot(x,Data[1][:,0,0],'--')
plt.plot(x,Data[1][:,1,1],'--')
plt.plot(x,Data[1][:,2,2],'--')
plt.plot(x,-Data[1][:,0,2],'--')
plt.xscale('log')
plt.yscale('log')
plt.show()
out = np.vstack([x, y, Data[1][:,0,0], Data[1][:,1,1], Data[1][:,2,2], Data[1][:,0,2]]).T
np.savetxt(output_folder+case+'.csv', out, newline='\n', delimiter=',', fmt='%10.5f')

