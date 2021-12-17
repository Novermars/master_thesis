
import torch
import matplotlib.pyplot as plt
from pylab import *
import pickle


from source.Calibration import CalibrationProblem
from source.DataGenerator import OnePointSpectraDataGenerator

####################################
### Configuration
####################################

case = 'tauNet_Kaimal_pen_reg' ### 'tauNet_Kaimal_noreg', 'tauNet_Kaimal_pen_reg', 'Mann_Kaimal', 'tauNet_Kaimal', 'tauNet_Kaimal_cross', 'tauNet_SimuYeo'
input_folder  = '/path/to/code/data/'
output_folder = '/path/to/figures/Convergence/'


####################################
### Import results
####################################
with open(input_folder+case+'.pkl', 'rb') as file:
    config, parameters, Data, _, loss_history = pickle.load(file)


####################################
### Convergence plot
####################################
x = np.arange(len(loss_history))
y = loss_history
plt.plot(x,y, 'o--')
# plt.xscale('log')
plt.yscale('log')
plt.show()
out = np.vstack([x, y]).T
np.savetxt(output_folder+'loss_convergence.csv', out, newline='\n', delimiter=',', fmt='%10.5f')

