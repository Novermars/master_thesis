
import torch
import matplotlib.pyplot as plt
from pylab import *
import pickle


from source.Calibration import CalibrationProblem
from source.DataGenerator import OnePointSpectraDataGenerator
from source.common import MannEddyLifetime

####################################
### Configuration
####################################

case = 'tauNet_Kaimal_pen_reg' ### 'tauNet_Kaimal_noreg', 'tauNet_Kaimal_pen_reg', 'Mann_Kaimal', 'tauNet_Kaimal', 'tauNet_Kaimal_cross', 'tauNet_SimuYeo'
input_folder  = '/path/to/code/data/' ## NOTE: CHANGE PATH
output_folder = '/path/to/figures/EddyLifetime_plot/' ## NOTE: CHANGE PATH


####################################
### Import results
####################################
with open(input_folder+case+'.pkl', 'rb') as file:
    config, parameters, Data, _, _ = pickle.load(file)

# config['learn_nu'] = False

pb = CalibrationProblem(**config)
pb.parameters = parameters


k_gd = torch.logspace(-3,3,50, dtype=torch.float64)
k_1  = torch.stack([k_gd, 0*k_gd, 0*k_gd], dim=-1)
k_2  = torch.stack([0*k_gd, k_gd, 0*k_gd], dim=-1)
k_3  = torch.stack([0*k_gd, 0*k_gd, k_gd], dim=-1)
k_4  = torch.stack([k_gd, k_gd, k_gd], dim=-1) / 3**(1/2)


####################################
### tau plot
####################################
tau_model1 = pb.OPS.EddyLifetime(k_1).detach().numpy()
tau_model2 = pb.OPS.EddyLifetime(k_2).detach().numpy()
tau_model3 = pb.OPS.EddyLifetime(k_3).detach().numpy()
tau_model4 = pb.OPS.EddyLifetime(k_4).detach().numpy()
tau_ref    = 3.9*MannEddyLifetime(0.59 * k_gd).detach().numpy()
plt.plot(k_gd, tau_model1, '-', label=r'$\tau_{model}(k_1)$')
plt.plot(k_gd, tau_model2, '-', label=r'$\tau_{model}(k_2)$')
plt.plot(k_gd, tau_model3, '-', label=r'$\tau_{model}(k_3)$')
plt.plot(k_gd, tau_model4, '-', label=r'$\tau_{model}(k,k,k)$')
plt.plot(k_gd, tau_ref,  '--', label=r'$\tau_{ref}=$Mann')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$k$')
plt.ylabel(r'$\tau$')
plt.grid(which='both')
plt.title('Eddy liftime')
plt.show()

out = np.vstack([k_gd, tau_ref, tau_model1, tau_model2, tau_model3, tau_model4 ]).T
np.savetxt(output_folder+'tau_plot.csv', out, newline='\n', delimiter=',', fmt='%10.5f')

