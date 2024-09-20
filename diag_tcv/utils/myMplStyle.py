#default
import matplotlib as mpl
from cycler import cycler

### Fontsize, figsize
mpl.rcParams["font.size"] = 16
mpl.rcParams["figure.figsize"] = [8,5]

### Color of lines and marker cyclers
# ~ mpl.rcParams["axes.prop_cycle"] = cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])).
# mpl.rcParams["axes.prop_cycle"] = (cycler('color', ['teal', 'firebrick', 'green', 'blue', 'indigo', 'brown', 'crimson', '#7f7f7f', '#bcbd22', '#17becf', 'red'])+
#                                     cycler('marker', ['o', '^', 's', 'D', 'X', '*', 'H', '<', '8', 'v', 'P']) )


mpl.rcParams["axes.prop_cycle"] = (cycler('color', ['blue', 'red', 'green','teal', 'firebrick', 'indigo', 'brown', 'crimson', '#7f7f7f', '#bcbd22', '#17becf']) )

# ~ mpl.rcParams['grid.linestyle'] = 'dotted'
mpl.rcParams['axes.grid'] = False

### plt.plot parameter
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams["lines.markersize"] = 8
mpl.rcParams["lines.markeredgecolor"] = 'black'
mpl.rcParams["lines.markeredgewidth"] = 1

### ticks parameters
mpl.rcParams['xtick.top'] = True
mpl.rcParams["xtick.bottom"] = True
mpl.rcParams["xtick.major.width"] = 2
mpl.rcParams["xtick.minor.width"] = 2
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["xtick.minor.size"] = 3
mpl.rcParams["xtick.color"] = "black"
mpl.rcParams['xtick.direction'] = 'in'


mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.left'] = True
mpl.rcParams["ytick.major.width"] = 2
mpl.rcParams["ytick.minor.width"] = 2
mpl.rcParams["ytick.major.size"] = 6
mpl.rcParams["ytick.minor.size"] = 3
mpl.rcParams["ytick.color"] = "black"
mpl.rcParams['ytick.direction'] = 'in'
