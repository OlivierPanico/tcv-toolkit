#%%
import scipy.io 

# data = scipy.io.loadmat('/home/panico/NoTivoli/astra_data/ASTRA_TCV_80376_0.400_1.935.mat')
data = scipy.io.loadmat('/home/panico/NoTivoli/astra_data/ASTRA_TCV_80328_0.400_2.000.mat')



# %%
from dataAnalysis.utils.plot_utils import plot_1d, plot_2d
import matplotlib.pyplot as plt
import numpy as np

rho = data['out']['RHOPSI'][0,0]
T = data['out']['T'][0,0][0]
QIEFF = data['out']['QIEFF'][0,0]
QEEFF = data['out']['QEEFF'][0,0]


print(rho.shape, T.shape, QIEFF.shape, QEEFF.shape)

# X,Y = np.meshgrid(rho[:,100], T)
# Z = np.transpose((QEEFF/(QEEFF+QIEFF)))
# fig, ax,_ = plot_2d(Z,X,Y, cmap=plt.cm.jet)
# ax.contour(X,Y,Z, [0.5], colors='white')

itime = 150
fig, ax = plot_1d([],[], grid=True)
# ax.plot(rho[:,itime], QEEFF[:,itime]/(QEEFF[:,itime]+QIEFF[:,itime]), label='QEEFF/(QEEFF+QIEFF)')
ax.plot(rho[:,itime], QEEFF[:,itime], label='QEEFF')
ax.plot(rho[:,itime], QIEFF[:,itime], label='QIEFF')
ax.set_title(f'Time: {T[itime]}')
ax.legend()


# %%
