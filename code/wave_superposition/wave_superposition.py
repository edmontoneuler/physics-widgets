import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib.widgets import Slider

#Position grid
x_min = -10
x_max = 10
dx = 0.01
x_vec = np.arange(x_min, x_max, dx)
ones = np.ones_like(x_vec)
#Time grid
t_min = 0
t_max = 10
t_vec = np.linspace(t_min, t_max, 10000)
dt = t_vec[1] - t_vec[0]
#Wave velocities
v1 = 1
v2 = -1
#Initial Parameters
x1 = -8
x2 = 5
w1 = 6
w2 = 4
h1 = 0.5
h2 = -0.5

def heaviside(arg):
    """Heaviside step function"""
    if arg >= 0:
        return 1
    else:
        return 0

def triangle(x, x0 = -8, w = 6, h = 0.5):
    """ Returns array defining right triangle waveform"""
    output = [ (heaviside(k - x0) - heaviside(k-(x0+w)))*h*(1-(k - x0)/w) for k in x]
    return np.array(output)

def box(x, x0 = 5, w = 4, h = -0.5):
    output = [ h*(heaviside(k - x0) - heaviside(k-(x0+w))) for k in x]
    return np.array(output)

#Initialization
fig, axs = plt.subplots(2,1)

xgrid_spacing = 1
ygrid_spacing = 0.1
axs[0].set_xticks(np.arange(x_min, x_max+xgrid_spacing, xgrid_spacing))
axs[0].set_yticks(np.arange(min(0, h1, h2), max(0, h1, h2)+ygrid_spacing, ygrid_spacing))
axs[1].set_xticks(np.arange(x_min, x_max+xgrid_spacing, xgrid_spacing))
axs[1].set_yticks(np.arange(min(0, h1, h2), max(0, h1, h2)+ygrid_spacing, ygrid_spacing))

axs[0].grid()
axs[1].grid()


wave1, = axs[0].plot(x_vec, triangle(x_vec), 'b', lw=2)
wave2, = axs[0].plot(x_vec, box(x_vec), 'r', lw=2)
pointA1, = axs[0].plot(x1, 0, 'ko')
pointA2, = axs[0].plot(x1, h1, 'ko')
pointA3, = axs[0].plot(x1+w1, 0, 'ko')
pointB1, = axs[0].plot(x2, 0, 'ko')
pointB2, = axs[0].plot(x2, h2, 'ko')
pointB3, = axs[0].plot(x2+w2, 0, 'ko')
pointB4, = axs[0].plot(x2+w2, h2, 'ko')
superwave, = axs[0].plot(x_vec, box(x_vec)+triangle(x_vec), 'k--', lw=1)
super_separate, = axs[1].plot(x_vec, box(x_vec)+triangle(x_vec), 'k', lw=3)

#Sliders
t_loc = plt.axes([0.25, 0.00, 0.50, 0.02])
t_amp = Slider(t_loc, 'Time (s)', t_min, t_max, valinit = t_min)

#Labels
axs[0].set_xlabel('x (m)')
axs[0].set_ylabel('y (m)')
axs[1].set_xlabel('x (m)')
axs[1].set_ylabel('y (m)')
axs[1].set_title('Wave Superposition')
axs[0].legend(['Wave 1', 'Wave 2', 'Superposition'])

def update(val):
    t = t_amp.val
    t_index = int((t-t_min)/dt)
    t = t_vec[t_index]
    
    W1 = triangle(x_vec-v1*t*ones)
    W2 = box(x_vec-v2*t*ones)
    wave1.set_ydata(W1)
    wave2.set_ydata(W2)
    superwave.set_ydata(W1+W2)
    super_separate.set_ydata(W1+W2)

    pointA1.set_data(x1+v1*t, 0)
    pointA2.set_data(x1+v1*t, h1)
    pointA3.set_data(x1+w1+v1*t, 0)
    pointB1.set_data(x2+v2*t, 0)
    pointB2.set_data(x2+v2*t, h2)
    pointB3.set_data(x2+w2+v2*t, 0)
    pointB4.set_data(x2+w2+v2*t, h2)

fig.canvas.draw_idle()
t_amp.on_changed(update)

plt.show()
