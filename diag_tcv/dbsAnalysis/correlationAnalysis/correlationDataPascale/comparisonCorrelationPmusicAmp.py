#%%
from diag_tcv.dbsAnalysis.correlationAnalysis.correlationAnalysis import CorrelationAnalysis
from dataAnalysis.utils.plot_utils import plot_1d, my_text, my_legend
#%%
a = CorrelationAnalysis(81069)
dircpmus = a.prepare_coherence_for_plot(10,x_norm='rho_s', load_from_pascale=True, retdata=True)
dircamp = a.prepare_coherence_for_plot(10,x_norm='rho_s', load_from_pascale=False, retdata=True)

#%%
fig, ax = plot_1d([],[], grid=True)
ax.plot(dircpmus['plateau0']['delta'], dircpmus['plateau0']['maxcorr'], marker='o', linestyle='None', label='pmusic')
ax.plot(dircamp['plateau0']['delta'], dircamp['plateau0']['maxcorr'], marker='^', linestyle='None', label='amp')
ax.set_yscale('log')
ax.set_ylim(0.05,1.1)

fig, ax = plot_1d([],[], grid=True)
ax.plot(dircpmus['plateau1']['delta'], dircpmus['plateau1']['maxcorr'], marker='o', linestyle='None', label='pmusic')
ax.plot(dircamp['plateau1']['delta'], dircamp['plateau1']['maxcorr'], marker='^', linestyle='None', label='amp')
ax.set_yscale('log')
ax.set_ylim(0.05,1.1)

fig, ax = plot_1d([],[], grid=True)
ax.plot(dircpmus['plateau2']['delta'], dircpmus['plateau2']['maxcorr'], marker='o', linestyle='None', label='pmusic')
ax.plot(dircamp['plateau2']['delta'], dircamp['plateau2']['maxcorr'], marker='^', linestyle='None', label='amp')
ax.set_yscale('log')
ax.set_ylim(0.05,1.1)
# %%
