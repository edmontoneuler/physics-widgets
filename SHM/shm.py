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

A = 1      #Amplitude
frequency = 1
w = 2*np.pi*frequency
wavelength = 1
v = wavelength*frequency
point_loc = 0

#Initialization
fig, axs = plt.subplots(3,1)

xgrid_spacing = 1
ygrid_spacing = 0.2
axs[0].set_xticks(np.arange(x_min, x_max+xgrid_spacing, xgrid_spacing))
axs[0].set_yticks(np.arange(-A, A+ygrid_spacing, ygrid_spacing))
axs[0].grid()
axs[1].set_xticks(np.arange(x_min, x_max+xgrid_spacing, xgrid_spacing))
axs[1].set_yticks(np.arange(-A*w, w*A+ygrid_spacing, w*ygrid_spacing))
axs[1].grid()
axs[2].set_xticks(np.arange(x_min, x_max+xgrid_spacing, xgrid_spacing))
axs[2].set_yticks(np.arange(-A*w*w, w*w*A+ygrid_spacing, w*w*ygrid_spacing))
axs[2].grid()


def position(x_vec):
    return np.sin(x_vec)

def velocity(x_vec, omega=w):
    return omega*np.cos(x_vec)

def acceleration(x_vec, omega=w):
    return -omega*omega*np.sin(x_vec)


position_wave, = axs[0].plot(x_vec, position(x_vec),'r', lw=2)
velocity_wave, = axs[1].plot(x_vec, velocity(x_vec), 'b', lw=2)
acceleration_wave, = axs[2].plot(x_vec, acceleration(x_vec), 'g', lw=2)
point_p, = axs[0].plot(point_loc, position(point_loc), 'ko')
point_v, = axs[1].plot(point_loc, velocity(point_loc), 'ko')
point_a, = axs[2].plot(point_loc, acceleration(point_loc), 'ko')

#Sliders
t_loc = plt.axes([0.25, 0.00, 0.50, 0.02])
t_amp = Slider(t_loc, 'Time (s)', t_min, t_max, valinit = t_min)

#Labels
#axs[0].set_xlabel('x (m)')
axs[0].set_ylabel('y (m)')
axs[0].set_title('Displacement')
#axs[1].set_xlabel('x (m)')
axs[1].set_ylabel('v (m/s)')
axs[1].set_title('Velocity')
axs[2].set_xlabel('x (m)')
axs[2].set_ylabel('a(m/s^2)')
axs[2].set_title('Acceleration')

def update(val):
    t = t_amp.val
    t_index = int((t-t_min)/dt)
    t = t_vec[t_index]

    position_wave.set_ydata(position(x_vec - v*t*ones))
    velocity_wave.set_ydata(velocity(x_vec - v*t*ones))
    acceleration_wave.set_ydata(acceleration(x_vec-v*t*ones))
    point_p.set_ydata(position(point_loc - v*t))
    point_v.set_ydata(velocity(point_loc - v*t))
    point_a.set_ydata(acceleration(point_loc - v*t))

    
    

fig.canvas.draw_idle()
t_amp.on_changed(update)

plt.show()
