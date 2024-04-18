#%%
from tcv_diag import TCVShot
import numpy as np
import matplotlib.pyplot as plt

def plot(rho, time, data, fmt='-', ax=None, twindow=[0, 2.5], rhowindow=[0.0, 1.4], yerr=None, errorband=False, errorband_alpha=0.2, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    ind = (time >= twindow[0]) & (time <= twindow[1])
    # time2d = np.tile(time, (rho.shape[1],1)).T
    # ind = (time2d >= twindow[0]) & (time2d < twindow[1]) & (rho >= rhowindow[0]) & (rho < rhowindow[1])

    t = time[ind]
    x = rho[ind]
    y = data[ind]

    if yerr is not None:

        yerr = yerr[ind]

        for i in range(x.shape[0]):
            
            col = f'C{i}'
            label = f'{t[i]:.2f} s'
            if errorband:
                l = ax.plot(x[i], y[i], fmt, label=label, color=col, **kwargs)
                # draw error bands:
                ax.fill_between(x[i], y[i]-yerr[i], y[i]+yerr[i], color=l[0].get_color(), alpha=errorband_alpha)
            else:
                l = ax.errorbar(x[i], y[i], yerr[i], fmt=fmt, label=label, color=col, **kwargs)

    else:
        
        for i in range(x.shape[0]):
            label = f'{t[i]:.2f} s'
            l = ax.plot(x[i], y[i], fmt, label=label, color=f'C{i}', **kwargs)
        
    return l

class CXRS(TCVShot):
    
    def __init__(self, shot, twindow=[0, 2.5]):
        super().__init__(shot)
        self.get_crxs()
        self.get_crxs_nofit()
        self.get_crxs_sys()
        
        self.shot = shot
        self.twindow = twindow
        
        
    def get_crxs_sys(self):
        """Supercharging the base class method"""
        
        for nsys in [1,2,3,4]:
            for quantity in ['vi', 'ni', 'ti']:
                
                node_path = r'\results::cxrs:acq_00' + f'{nsys}:{quantity.upper()}'
                node_path_err = r'\results::cxrs:acq_00' f'{nsys}:{quantity.upper()}:ERR'
                key = f'{quantity}_sys{nsys}'
                try:
                    val = self.tree.getNode(node_path)
                    val_err = self.tree.getNode(node_path_err)
                    
                    setattr(self, f'{key}_node', val)
                    setattr(self, f'{key}', val.data())
                    # setattr(self, f'{key}_err', np.nan * np.ones_like(val.data()))
                    setattr(self, f'{key}_err', val_err.data())
                    setattr(self, f'time_{key}', val.getDimensionAt(1).data())
                    setattr(self, f'rho_{key}', val.getDimensionAt(0).data())
                    
                except Exception as e:
                    print(e)
                    print(f'No "{quantity}_sys{nsys}" (on node path "f{node_path}") node for shot {self.shot}')
        # try:               
        #     super().get_crxs_sys()
        # except Exception as e:
        #     print(e)
        #     print(f'No "vtor_sys{nsys}" (on node path "f{node_path}") node for shot {self.shot}')
  
    def plot(self,  fmt='-', quantity='vtor_raw', twindow=None, ax=None, **plot_kwargs):
        
        unit_f = 1
        
        if twindow is None:
            twindow = self.twindow
            
        if quantity == 'vtor_raw':
            rho = self.rho_vtor_raw
            time = self.time_vtor_raw
            y = self.vtor_raw
            yerr = self.vtor_err_raw
            ylabel = 'Vtor [km/s]'
        elif quantity == 'vpol_raw':
            rho = self.rho_vpol_raw
            time = self.time_vpol_raw
            y = self.vpol_raw
            yerr = self.vpol_err_raw
            rho = self.rho_c3
            time = self.time_c3
            y = self.vtor_c3
            yerr = self.vpol_err_raw
            # yerr = self.vpol_err_raw
            
            ylabel = 'Vpol [km/s]'
            
        elif '_sys' in quantity:
            nsys = int(quantity[-1])
            _quantity = quantity.split('_')[0] # one of ['vi', 'ni', 'ti']
            
            try:
                        
                rho = getattr(self, f'rho_{_quantity}_sys{nsys}')
                time = getattr(self, f'time_{_quantity}_sys{nsys}')
                y = getattr(self, f'{_quantity}_sys{nsys}')
                yerr = getattr(self, f'{_quantity}_sys{nsys}_err')
                
                ylabel = f'{_quantity} sys{nsys} [km/s]' if _quantity == 'vi' else f'{_quantity} sys{nsys} [1e19 m^-3]' if _quantity == 'ni' else f'{_quantity} sys{nsys} [eV]'
                
                if _quantity == 'ni':
                    unit_f = 1e-19
                    
            except AttributeError as e:
                print(e)
                return None, None
        
        elif quantity == 'vtor_proffit':    
            
            rho = self.rho_vtor
            time = self.time_vtor
            y = self.vtor
            yerr = self.vtor_err
            ylabel = 'Vtor (proffit) [km/s]'
            
        elif quantity == 'vpol_proffit':
            
            rho = self.rho_vpol
            time = self.time_vpol
            y = self.vpol
            yerr = self.vpol_err
            ylabel = 'Vpol (proffit) [km/s]'
        
        elif quantity == 'ni_proffit':
            
            rho = self.rho_ni
            time = self.time_ni
            y = self.ni
            yerr = self.ni_err
            ylabel = 'Ni (proffit) [1e19 m^-3]'
            unit_f = 1e-19
            
        elif quantity == 'ti_proffit':
            
            rho = self.rho_ti
            time = self.time_ti
            y = self.ti
            yerr = self.ti_err
            
            ylabel = 'Ti (proffit) [eV]'
            
        else:
            raise NotImplementedError(f'quantity {quantity} not implemented yet.')
        
        twindow[0] = max(twindow[0], time[0])
        twindow[1] = min(twindow[1], time[-1])
     
        l = plot(rho, time, y * unit_f, fmt, yerr=yerr * unit_f, ax=ax, twindow=twindow, **plot_kwargs)
        
        ax.set_ylabel(ylabel)
        ax.set_xlabel(r'$\rho_\psi$')
        
        data = {'rho': rho, 'time': time, 'y': y, 'yerr': yerr}
        
        return l, data


def cxrs_overview(shot, twindow=None, axs=None, cxrs_obj=None, show_sys=[1,2,3,4], show_legend=True, show_proffit=True, make_title=True, **plot_kwargs):

    
    sys_markers = {
        1: 'x',
        2: 'o',
        3: 's',
        4: 'd'
    }
    
    quantity_dict = {
        'vi_sys_tor': {'ylabel': r'$v_\varphi$ [km/s]', # (sys 1 + 2)
                       'nsys': [1,2]},
        'vi_sys_pol': {'ylabel': r'$v_\theta$ [km/s]', # (sys 3 + 4)
                          'nsys': [3,4]},
        'ni_sys': {'ylabel': r'$n_C$ [1e19 m$^-3$] ',# (sys 1-3)
                     'nsys': [1,2,3,4]},
        'ti_sys': {'ylabel': r'$T_C$ [eV] ', # (sys 1-4)
                    'nsys': [1,2,3,4]}
    }
    
    from matlabtools import Struct
    plot_items = Struct()
    legends1 = []
    legends2 = []
    
    errorband = plot_kwargs.pop('errorband', True)
    errorband_alpha = plot_kwargs.pop('errorband_alpha', 0.2)
    default_plot_kwargs = dict(alpha=0.5, ms=4, capsize=0.2, lw=1)
    
    if twindow is None:
        twindow = [0, 2.5]
    
    if cxrs_obj is None:
        cxrs= CXRS(shot, twindow=twindow)
    else:
        cxrs = cxrs_obj
        
    if axs is None:
        fig, axs = plt.subplots(2, 2, sharex=True, sharey=False, figsize=(8,6))

    
    # proffit curves:
    if show_proffit:
        for i, (quantity, ax) in enumerate(zip(['vtor', 'vpol', 'ni', 'ti'],axs.flat)):
            
            l, data = cxrs.plot('-', ax=ax, quantity=f'{quantity}_proffit', errorband=errorband, errorband_alpha=errorband_alpha)
            ax.grid(True)
            
            if show_legend:
                leg1 = ax.legend(loc='upper right', handlelength=0.5, fontsize=8, ncol=2) # keep the reference to the first legend  
                legends1.append(leg1)
            
            plot_items[f'{quantity}_proffit'] = {'plot_item': l, 'data': data}
    
    # raw scatter plots from individual CXRS systems
    for i, (quantity, ax) in enumerate(zip(['vi_sys_tor', 'vi_sys_pol', 'ni_sys', 'ti_sys'],axs.flat)):

        for key, val in default_plot_kwargs.items():
            plot_kwargs[key] = plot_kwargs.get(key, val)
            
        _nsys = quantity_dict[quantity]['nsys']
        _ylab = quantity_dict[quantity]['ylabel']

        data = None # placeholder to store the data from the different systems
        for nsys in _nsys:
            if nsys not in show_sys:
                continue
            
            l, data = cxrs.plot(sys_markers[nsys], ax=ax, quantity=f'{quantity}_sys{nsys}', **plot_kwargs)
            
            # if data is None:
            #     data = _data
            # else:
            #     data['rho'] = np.concatenate([data['rho'], _data['rho']])
            #     # data['time'] = np.concatenate([data['time'], _data['time']])
            #     data['y'] = np.concatenate([data['y'], _data['y']])
            #     data['yerr'] = np.concatenate([data['yerr'], _data['yerr']])
                
  
        ax.set_ylabel(_ylab)
        plot_items[f'{quantity}_sys{nsys}'] = {'plot_item': l, 'data': data, 'ax':ax}
        
            # for nsys in [1,2]:
            #     if nsys not in show_sys:
            #         continue
            #     l, data = cxrs.plot(sys_markers[nsys], ax=ax, quantity=f'vi_sys{nsys}', **plot_kwargs)
                    
            # ax.set_ylabel(r'$v_\varphi$ [km/s] (sys 1 + 2)')
            
        # elif quantity == 'vi_sys_pol':
        #     for nsys in [3,4]:
        #         if nsys not in show_sys:
        #             continue
        #         l, data = cxrs.plot(sys_markers[nsys], ax=ax, quantity=f'vi_sys{nsys}', **plot_kwargs)
        #     ax.set_ylabel(r'$v_\theta$ [km/s] (sys 3 + 4)')
                
        # elif quantity == 'ni_sys':
        #     for nsys in [1,2,3,4]:
        #         if nsys not in show_sys:
        #             continue
                
        #         l, data = cxrs.plot(sys_markers[nsys], ax=ax, quantity=f'ni_sys{nsys}', **plot_kwargs)
        #     ax.set_ylabel(r'$n_i$ [1e19 m^-3] (sys 1-3)')
        # elif quantity == 'ti_sys':
        #     for nsys in [1,2,3,4]:
                
        #         if nsys not in show_sys:
        #             continue
                
        #         l, data = cxrs.plot(sys_markers[nsys], ax=ax, quantity=f'ti_sys{nsys}', **plot_kwargs)
        #     ax.set_ylabel(r'$T_i$ [eV] (sys 1-4)')
        
        # else:
        #     l, data = cxrs.plot('x', ax=ax, quantity=quantity, **plot_kwargs)
            
    
        # add a nother custom legend:
        from matplotlib.lines import Line2D
        handles = []
        for nsys, marker in sys_markers.items():
            if nsys not in show_sys or nsys not in _nsys:
                continue
            l, = ax.plot([], [], marker, color='k', label=f'sys {nsys}')
            handles.append(l)
            
        leg2 = ax.legend(handles=handles, handlelength=0.5, fontsize=8, framealpha=0.6, loc='upper left')
        legends2.append(leg2)

    if show_legend:
        for ax,leg in zip(axs.flat, legends1):
            ax.add_artist(leg)
            
    if make_title:
        ax.set_title(f'TCV #{cxrs.shot}, t={twindow[0]:.2f}-{twindow[1]:.2f} s')
    
    ax.set_xlim(left=0.6, right=1.15)
    plt.tight_layout()
    
    plot_items['legends_times'] = legends1
    plot_items['legends_nsys']  = legends2
    
    return plot_items

def check_vtor_vpol_consistency(shot=79349, twindow=[1.0,1.2], cxrs_obj=None):
    """Just to check for consistency (note the wiki page used to be erroneous - it said system 2 where it should be 3 for vpol_raw)
    cxrs= TCVCXRS(shot, twindow=twindow)"""
    if cxrs_obj is None:
        cxrs= CXRS(shot, twindow=twindow)
    else:
        cxrs = cxrs_obj
    
    fig,axs = plt.subplots(1,2, sharex=True, sharey=False)
    for i, (quantity, ax) in enumerate(zip(['vtor_raw', 'vpol_raw'],axs)):
            l = cxrs.plot('x', ax=ax, quantity=quantity, alpha=0.5)
    for i, (quantity, ax) in enumerate(zip(['vi_sys1', 'vi_sys3'],axs)):
            l = cxrs.plot('o', ax=ax, quantity=quantity, alpha=0.5)

#%%
if __name__ == '__main__':
    cxrs_obj = CXRS(79335, twindow=[0.5,1.0])
    #%%
    plt_items = cxrs_overview(cxrs_obj.shot, cxrs_obj=cxrs_obj)
        
    #%%
    check_vtor_vpol_consistency()

#%%



cxrs_obj = CXRS(80257)#, twindow=[0,2.5])
    
rho_ti = cxrs_obj['rho_ti']
time_ti = cxrs_obj['time_ti']
ti = cxrs_obj['ti']
ti_err = cxrs_obj['ti_err'] 

#%%
#print("rho ti position = " , rho_ti[0, 88])
fig, ax = plt.subplots()
a,b = cxrs_obj.plot(quantity='ti_proffit', ax=ax,rhowindow=[0.5, 1], twindow=[0.9,1], errorband=True)
plt.xlim(0.2,1)
plt.legend()

#%%
plt.plot(time_ti, ti[:,0], marker='+')

#%%

plt.figure()
plt.pcolormesh(rho_ti[10,:], time_ti, ti, cmap=plt.cm.jet)


