## REMEMBER TO ACTIVATE `GPU_optimization (aGPUo)`
#import tensorflow
import torch
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scienceplots
import imageio.v2 as imageio
import glob

from matplotlib import animation
from matplotlib.animation import PillowWriter

plt.style.use(['science', 'notebook'])
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# from tensorflow.python.client import device_lib
# print(device_lib.list_local_devices())


# Parameters | Remember it is the normalized analysis

L  = 100    # Domain of the solution 0 <= x <= L  (in Debye lengths)
N  = 25000  # Number of electrons
J  = 1000   # Number of grid-points
vb = 5      # Beam velocity
n0 = N/L    # ion number density
dx = L/J

dt = 0.1    # time step  (in inverse plasma frequencies)
t_max = 150  # such that 0 <= t <= t_max
timesteps = int(t_max / dt)


# Check input parameters make sence:
if (N < 1) | (J < 2) | (L <= 0.) | (vb <= 0.) | (dt <= 0.) | (t_max <= 0.) | ((int (t_max / dt) / 10) < 1):
    print("Error - invalid input parameters")
   
   
# Initial positions
r0 = np.random.random(N)*L


# Initial velocities: utilizing REJECTION SAMPLING
def sample_velocity(v_b, tails_factor):
    f_max = 0.5 * (1.0 + np.exp(-2.0 * vb**2))
    vmin, vmax = -tails_factor*v_b, tails_factor*v_b
    
    velocities = np.zeros(N)
    i = 0
    while i < N:
        v_ = np.random.uniform(vmin, vmax) 
        beam_shape = 0.5 * (np.exp(-0.5 * (v_ - v_b)**2) + np.exp(-0.5 * (v_ + v_b)**2))
        gamma = np.random.rand()*f_max
        if gamma <= beam_shape:
            velocities[i] = v_
            i += 1
    return velocities

v0 = sample_velocity(vb, 4) 
velocity_tags = abs(v0 - vb) < abs(v0 + vb) # False = left going


def compute_density(r_):
    ne_ = np.zeros(J)
    j_ = (r_/dx).astype(int)         # size N
    fractional_offset = (r_/dx) - j_ # distance from left grid point (normalized)
    
    #for i in range(N):
    #    left,right = j_[i]%J, (j_[i]+1)%J               # particle's indexes w.r.t. j
    #    ne_[right] += fractional_offset[i] / dx          # weight to left
    #    ne_[left ] += (1 - fractional_offset[i]) / dx    # weight to right
    
    # REMOVING LOOP TO IMPLEMENT
    np.add.at(ne_, j_%J, (1 - fractional_offset) / dx)
    np.add.at(ne_, (j_+1)%J, fractional_offset / dx)

    
def poisson_solver(ne_):
    # 1. Compute normalized charge density ρ(x) = ne/n0 - 1
    rho = ne_ / n0 - 1

    # 2. Compute the FFT of ρ(x)
    # This returns complex coefficients for Fourier modes:
    # rho_hat[0] → k=0 (mean)
    # rho_hat[1:J//2] → positive frequencies
    # rho_hat[J//2+1:] → negative frequencies (wrapped around)
    rho_hat = np.fft.fft(rho)

    # 3. Define wavenumbers k = 2πn / L for n = 0 to J-1
    # These are scaled correctly to match the physical domain
    k = 2 * np.pi * np.fft.fftfreq(J, d=L/J)

    # 4. Allocate array for φ̂(k) — the Fourier coefficients of φ
    phi_hat = np.zeros(J, dtype=complex)

    # 5. Solve φ̂(k) = -ρ̂(k) / (k^2), avoiding k=0 to prevent division by zero
    # This corresponds to solving: d²φ/dx² = ρ → φ̂ = -ρ̂ / k²
    # We set φ̂[0] = 0 to ensure φ has no DC (mean) component
    phi_hat[1:] = -rho_hat[1:] / (k[1:]**2)
    phi_hat[0] = 0.0  # No DC offset

    # 6. Inverse FFT to get back to φ(x) in real space
    # The real part is taken because the physical φ(x) must be real
    phi = np.fft.ifft(phi_hat).real

    return phi


def electric_field(phi_):
    E = np.zeros(J)
    for j in range(J):
        E[j] = ( phi_[(j-1)%J] - phi_[(j+1)%J])/(2*dx)
    return E


def particle_electric_field(r_, E_):
    E_particles = np.zeros(N)
    j_ = (r_/dx).astype(int)  # size N
    y_ = (r_/dx) - j_         # fractional offset
    
    # ACTUALLY, A LOOP IS NOT NEEDED HERE, IT CAN BE IMPROVED WITH NUMPY'S VECTORIZATION
    for i in range(N):
        left,right = j_[i]%J, (j_[i]+1)%J   # particle's indexes w.r.t. j
        E_particles[i] = (1-y_[i])*E_[left] + y_[i]*E_[right]  # Linear interpolation 
    return E_particles


def derivative_r_v(y_):
    """
    Returns dy/dt, with f = [r0,r1, ..., rN,v0,v1,...,vN]
    Such that dy/dt = [dr/dt, dv/dt] = [v(r_i), -E(r_i)]
    For notation dy/dt --> dy [tuple]
    """
    
    r_ = y_[:N]
    v_ = y_[N:] 
    
    # 1. Compute the density
    ne = compute_density(r_)
    
    # 2. Solve Poisson's equation
    phi = poisson_solver(ne)
    
    # 3. Compute the Electric field in the grid
    E_grid = electric_field(phi)
    
    # 4. Compute the Electric field per particle
    E_particles = particle_electric_field(r_,E_grid)
    
    return np.concatenate([v_, -E_particles])
    
    
def RungeKutta4(y_):
    """
    dy = dy/dt = [dr/dt, dv/dt] = [v(r_i), -E(r_i)]
    """
    
    # f(x,y)    
    k1 = derivative_r_v( y_)
    k2 = derivative_r_v( y_ + 0.5*dt*k1)
    k3 = derivative_r_v( y_ + 0.5*dt*k2)
    k4 = derivative_r_v( y_ + dt*k3)
    
    return y_ + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)


sim_params = f"L{L}N{N}J{J}vb{vb}dt{dt}tmax{t_max}"
location = "Images/frames/" + sim_params
os.makedirs(location, exist_ok=True)
def save_phase_space_plot(r, v, step):
    # Create output directory for creation of a gif
    plt.figure(figsize=(8, 4))
    plt.plot(r[velocity_tags], v[velocity_tags], '.', markersize=1, alpha=0.3, color='r')
    plt.plot(r[~velocity_tags], v[~velocity_tags], '.', markersize=1, alpha=0.3, color='k')
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    plt.title(f"Phase Space - time {round(step*dt,2)}")
    plt.xlim(0, L)
    plt.ylim(-2*vb, 2*vb)
    plt.grid(True, alpha=0.5)
    filename = location + f"/frame_{step:04d}.png"
    plt.savefig(filename, dpi=150)
    plt.close()
    
def create_gif(output_path=location+"simulation_.gif", fps=10):
    frames = sorted(glob.glob(location+"/frame_*.png"))
    images = [imageio.imread(f) for f in frames]
    imageio.mimsave(output_path, images, fps=fps)


y = np.concatenate([r0,v0])
for step in range(timesteps):
    y = RungeKutta4(y)
    # Enforce periodic boundary conditions
    y[:N] = y[:N] % L
    # Save simulation to create a gif
    save_phase_space_plot(y[:N], y[N:], step)
    
create_gif()
    
