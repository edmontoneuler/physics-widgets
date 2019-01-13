import numpy as np
import matplotlib.pyplot as plt
import scipy 
import scipy.io as spio
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib.widgets import Slider
import sys

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key], spio.matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
        return d

    def _todict(matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, spio.matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, spio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

material_prefix = sys.argv[1]
data = loadmat(material_prefix + '_plot_Ek.mat')
eposts = data['eposts']
nbins = data['nbins']

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cm = plt.get_cmap("jet")

conduction_loc = plt.axes([0.25, 0.05, 0.50, 0.02])
conduction_amp = Slider(conduction_loc, 'Upper bound (eV)', eposts[0], eposts[-1], valinit=0.20)
valence_loc = plt.axes([0.25, 0.00, 0.50, 0.02])
valence_amp = Slider(valence_loc, 'Lower bound (eV)', eposts[0], eposts[-1], valinit=-0.05)

def lower_bin_index(energy):
    for k in range(len(eposts)):
        if energy > eposts[k]:
            pass
        else:
            return k-1
             
def upper_bin_index(energy):
    for k in range(len(eposts)):
        if energy > eposts[k]:
            pass
        else:
            return k-1

def get_scatter_data (bin_list):
    X = []
    Y = []
    Z = []
    E = []


    for k in bin_list:
        temp_bin = data['ebins'][k].kpoints
        if np.shape(temp_bin)[0] != 0:
            n, _ = np.shape(temp_bin)
            for k in range(n):
                Y.append(temp_bin[k][0])
                Z.append(temp_bin[k][1])
                X.append(temp_bin[k][2])
                E.append(temp_bin[k][3])
    X = np.array(X)
    Y = np.array(Y)
    Z  = np.array(Z)
    E = np.array(E)
    return X, Y, Z, E
 
X0, Y0, Z0, E0 = get_scatter_data(np.arange(44, 48))
norm = plt.Normalize(min(E0), max(E0))
colors = cm(norm(E0))
scat = ax.scatter(X0, Y0, Z0, c=colors)

material_prefix = material_prefix[:-5]
def update(val):
    e_lower = valence_amp.val
    e_upper = conduction_amp.val
    bin_list = np.arange(lower_bin_index(e_lower), upper_bin_index(e_upper))
    X, Y, Z, E = get_scatter_data(bin_list)
    e_min, e_max = min(E), max(E)
    norm = plt.Normalize(e_min, e_max)
    colors = cm(norm(E))
    ax.clear()
    scat = ax.scatter(X, Y, Z, c=colors)
    
    ax.set_xlabel('k_x')
    ax.set_ylabel('k_y')
    ax.set_zlabel('k_z')
    ax.set_title( material_prefix + ' Bandstructure')

    fig.canvas.draw_idle()

valence_amp.on_changed(update)
conduction_amp.on_changed(update)
plt.show()
