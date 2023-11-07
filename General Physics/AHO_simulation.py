import numpy as np
import scienceplots
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Please, notice this code is strongly based on the free course 
# "Stochastic Processes: Data Analysis and Computer Simulation"
# found at EDX. 


plt.style.use(['science', 'notebook'])

sys_dim = 2  # system dimension (x,y)
steps   = 1500 # step number
R_ = np.zeros(sys_dim) # current position (x,y)
V_ = np.zeros(sys_dim) # current velocity (x,y)
Rt = np.zeros([sys_dim,steps]) # position for all times
Vt = np.zeros([sys_dim,steps]) # velocity for all times
E  = np.zeros(steps) # Energy for all times
t  = np.zeros(steps) # Time array

def init_animation():
    """
    Definition of all graphical elementes to be used
    """
    particles.set_data([],[])
    line.set_data([],[])
    title.set_text(r'')
    return particles, line, title

def animate_DHO(i):
    global R_, V_, Rt, Vt, R, t
    V_ = (1-ζ*Δt/m)*V_  -  k*(np.linalg.norm(R_)**2)*R_*Δt/m
    R_ = R_ + V_*Δt
    Rt[0:sys_dim, i] = R_
    Vt[0:sys_dim, i] = V_
    t[i] = i*Δt
    E[i] = 0.5*m*(np.linalg.norm(V_)**2) + 0.5*k*(np.linalg.norm(R_)**2)
    particles.set_data(R_[0],R_[1])                       # current position
    line.set_data(Rt[0,0:i], Rt[1,0:i])   # add latest position
    title.set_text(r'$t={0:.2f}, E_t = {1:.3f}$'.format(i*Δt, E[i]) )
    return particles, line, title

# particle's mass, spring, and frinction constant
m, k, ζ = 1.0, 1.0, 0.005
# initial condition
R_[0],R_[1] = 1.0, 1.0
V_[0],V_[1] = 1.5, 0.0
Δt = 0.1*np.sqrt(k/m)
# area for drawing
box_size = 5

# Setting up the figure for animation
fig, ax = plt.subplots(figsize=(7.5,7.5))
ax = plt.axes( xlim=(-box_size/2, box_size/2), ylim=(-box_size/2, box_size/2 ))
particles, = ax.plot([],[], 'ko', ms=10)
line,  = ax.plot([], [], lw=1, color="red")
title  = ax.text(0.5, 1.05, r'', transform=ax.transAxes, va='center' )
plt.grid(alpha=0.5)


animate_particle = animation.FuncAnimation(fig, animate_DHO, init_func=init_animation,
                                           frames=steps, interval=5, blit=True, repeat=False)


writergif = animation.PillowWriter(fps=30) 
animate_particle.save("AHO_zeta0_005.gif", writer=writergif)

print("hello")