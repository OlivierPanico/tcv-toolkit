#%%
import matplotlib.pyplot as plt
import numpy as np
from DBS.io.interface import DataInterface
from DBS.processing.sigprocessing import init_specobjs, show_spec, make_title, OutputWrapper
# a sweep that has been treated already:
machine, shot, channelval, isweep = 'tcv', 80376, 2, 18

data_interface = DataInterface(shot, isweep, channelval, machine=machine)
output_wrapper = OutputWrapper(data_interface) # automatically loads the existing data from file if any
ifreqs = data_interface._get_freq_choice('all')


# initialize the specobjs:
specobjs = init_specobjs(data_interface, ifreqs)

# check if some of the specobjs already have been treated (validated or rejected), in which case they are overwritten:
ifreqs_treated = output_wrapper.header.ifreqs_treated
if len(ifreqs_treated) > 0:
    print(f'Found {len(output_wrapper.header.ifreqs_treated)} specobjs already treated.')
        
    for i in ifreqs_treated:
        if i in ifreqs:
            k = (np.where(i==ifreqs))[0][0]
            print(k)
            specobjs[k] = output_wrapper.merged_specobj.specobjs[i-1]

output_wrapper.update_specobjs(specobjs)


for s in specobjs[:1]:
    # perform_specobj_fits(s)
    fig, ax = plt.subplots()
    plot_dict = show_spec(s, ax=ax, verbose=True)
    make_title(ax, data_interface, s.header.ifreq)
    
    
# %%
s=output_wrapper.specobjs[13]

fig, ax = plt.subplots()
from DBS.processing.fit_utils import fitfuncs

include_mask = s.include_mask

xscale = s.xscale
yscale = s.yscale

xdata = s.f / xscale
ydata = s.P / yscale
xdata_masked = s.f[include_mask] / xscale
ydata_masked = s.P[include_mask] / yscale

curve_type = 'taylor'
curve_func   = fitfuncs[curve_type]
curve_params = s.fit_params[curve_type]
xfit = np.linspace(np.min(xdata), np.max(xdata), len(xdata))#1000)
yfit = curve_func(xfit, *curve_params, dt=s.dt * xscale)

x = xfit * xscale
y = yfit * yscale + s.P_noise


l = ax.plot(xdata_masked * xscale / 1e6, 10 * np.log10(ydata_masked * yscale), label=f'raw')
l = ax.plot(x / 1e6, 10 * np.log10(y), label=f'{curve_type} fit')
ax.set_xlabel('f [MHz]')
ax.set_ylabel('PSD [dB]')
ax.legend()

# l = ax.plot(xdata_masked, 10 * np.log10(ydata_masked), label=f'raw')
# l = ax.plot(xfit, 10 * np.log10(yfit), label=f'{curve_type} fit')
# ax.set_xlabel('f [MHz]')
# ax.set_ylabel('PSD [dB]')
# ax.legend()


# %%
#Goal here is to take ydata and create a ydata_filled so that the excluded zones (include_mas=False) are filled with the fit curve
# ydata: list of y-values
# include_mas: list of boolean values
# fit_curve: function that calculates y-value from x-value

ydata_filled = []
for i in range(len(xdata)):
    if not include_mask[i]:
        ydata_filled.append(y[i])
    else:
        ydata_filled.append(ydata[i]*yscale)
ydata_filled=np.array(ydata_filled)

plt.figure()
plt.plot(s.f/xscale, 10*np.log10(s.P/yscale))
plt.plot(xdata_masked, 10*np.log10(ydata_masked))
plt.plot(xdata,10*np.log10(ydata_filled/yscale))

#%%
# spectrum_raw_ch1 = ydata
# spectrum_filled_ch1 = ydata_filled/yscale


spectrum_raw_ch2 = ydata
spectrum_filled_ch2 = ydata_filled/yscale


# %%
plt.figure()
plt.plot(xdata, spectrum_raw_ch1)
plt.plot(xdata, spectrum_raw_ch2)
plt.plot(xdata, spectrum_filled_ch1)
plt.plot(xdata, spectrum_filled_ch2)
plt.yscale('log')


#%%

psd_raw_ch1=spectrum_raw_ch1*np.conjugate(spectrum_raw_ch1)
psd_raw_ch2=spectrum_raw_ch2*np.conjugate(spectrum_raw_ch2)
psd_filled_ch1=spectrum_filled_ch1*np.conjugate(spectrum_filled_ch1)
psd_filled_ch2=spectrum_filled_ch2*np.conjugate(spectrum_filled_ch2)

cpsd_raw=spectrum_raw_ch1*np.conjugate(spectrum_raw_ch2)
cpsd_filled=spectrum_filled_ch1*np.conjugate(spectrum_filled_ch2)

coherence_raw=abs(cpsd_raw)**2/(psd_raw_ch1*psd_raw_ch2)
coherence_filled=abs(cpsd_filled)**2/(psd_filled_ch1*psd_filled_ch2)

plt.figure()
# plt.plot(xdata, psd_raw_ch1)
plt.plot(xdata, psd_filled_ch1)
plt.plot(xdata, psd_filled_ch2)
plt.yscale('log')


plt.figure()
plt.plot(xdata, np.abs(cpsd_raw))
plt.plot(xdata, np.abs(cpsd_filled))
# plt.yscale('log')

plt.figure()
plt.plot(xdata, coherence_raw)  
plt.plot(xdata, coherence_filled)

# %%
import numpy.fft as fft

def get_psd(cpsd, dt):
    return fft.ifft(cpsd).real/dt

# cpsd = cpsd_raw
dt = s.dt
psd_raw = get_psd(cpsd_raw, dt)
psd_filled = get_psd(cpsd_filled, dt)

plt.figure()
plt.plot(xdata, psd_raw)
plt.plot(xdata, psd_filled)
plt.yscale('log')

# %%
