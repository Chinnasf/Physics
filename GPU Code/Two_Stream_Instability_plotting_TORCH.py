import torch
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scienceplots
import imageio.v2 as imageio
import glob
import time

from matplotlib import animation
from matplotlib.animation import PillowWriter

plt.style.use(['science', 'notebook'])
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


torch.cuda.synchronize()
start_time = time.time()

L  = 100     # Domain of the solution 0 <= x <= L  (in Debye lengths)
N  = 250000  # Number of electrons
J  = 1000    # Number of grid-points
vb = 5       # Beam velocity
n0 = N/L     # ion number density
dx = L/J

dt = 0.1     # time step  (in inverse plasma frequencies)
t_max = 60  # such that 0 <= t <= t_max
timesteps = int(t_max / dt)

# FOR THIS SETTING:
# Elapsed: 404.86 s (t_max = 100)
# Elapsed: 252.56 s (t_max = 60)

# Check input parameters make sence:
if (N < 1) | (J < 2) | (L <= 0.) | (vb <= 0.) | (dt <= 0.) | (t_max <= 0.) | ((int (t_max / dt) / 10) < 1):
    print("Error - invalid input parameters")


def sample_velocity(v_b, tails_factor, batch_size=100000):
    f_max = 0.5 * (1.0 + np.exp(-2 * v_b**2))
    vmin, vmax = -tails_factor * v_b, tails_factor * v_b

    velocities = torch.empty(N, device=device)
    filled = 0

    while filled < N:
        v_ = (vmax - vmin) * torch.rand(batch_size, device=device) + vmin
        beam_shape = 0.5 * (torch.exp(-0.5 * (v_ - v_b)**2) + torch.exp(-0.5 * (v_ + v_b)**2))
        gamma = torch.rand(batch_size, device=device) * f_max
        accepted = v_[gamma <= beam_shape]

        num_to_fill = min(accepted.shape[0], N - filled)
        velocities[filled:filled+num_to_fill] = accepted[:num_to_fill]
        filled += num_to_fill

    return velocities


def compute_density(r_):
    ne_ = torch.zeros(J, device=device)
    j_  = torch.floor(r_ / dx).long()
    f   = (r_ / dx) - j_ # fractional offset
    
    # Right: weight = f / dx
    ne_.scatter_add_(0, (j_+1)%J, f / dx)
    
    # Left: weight = (1 - f) / dx
    ne_.scatter_add_(0, j_%J, (1 - f) / dx)
    
    return ne_



def poisson_solver(ne_):
    rho = ne_ / n0 - 1 # normalized charge density
    
    # FFT of rho (torch.fft returns complex tensor)
    rho_hat = torch.fft.fft(rho)
    
    # Wavenumbers (same as np.fft.fftfreq)
    k = 2 * np.pi * torch.fft.fftfreq(J, d=L/J).to(device)
    
    # Avoid division by zero
    phi_hat = torch.zeros_like(rho_hat)
    phi_hat[1:] = -rho_hat[1:] / (k[1:]**2)
    phi_hat[0] = 0.0 + 0.0j
    
    # Inverse FFT to get potential phi(x)
    phi = torch.fft.ifft(phi_hat).real  # result should be real-valued
    
    return phi


def electric_field(phi_):
    # Use periodic shifting
    phi_left  = torch.roll(phi_, shifts=1, dims=0)
    phi_right = torch.roll(phi_, shifts=-1, dims=0)
    E_grid = (phi_left - phi_right) / (2 * dx)
    return E_grid



def particle_electric_field(r_, E_):
    j_ = torch.floor(r_ / dx).long()
    f_ = (r_ / dx) - j_
    E_particles = (1 - f_)*E_[ j_ % J ] + f_*E_[ (j_ + 1) % J ]
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
    
    return torch.cat([v_, -E_particles])


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



# Initial positions
r0 = torch.rand(N, device=device) * L
v0 = sample_velocity(vb,4)
y  = torch.cat([r0, v0])

# False = left going | used for plotting
velocity_tags = (abs(v0 - vb) < abs(v0 + vb)).cpu().numpy()

sim_params = f"L{L}N{N}J{J}vb{vb}dt{dt}tmax{t_max}"
location = "Images/frames/" + sim_params
os.makedirs(location, exist_ok=True)

def save_phase_space_plot(r_, v_, step):
    # Create output directory for creation of a gif

    r = r_.cpu().numpy()
    v = v_.cpu().numpy()

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


#############
#
#  MAIN LOOP 
#
#############


for step in range(timesteps):
    y = RungeKutta4(y)
    # Enforce periodic boundary conditions
    y[:N] = y[:N] % L
    # Save simulation to create a gif
    save_phase_space_plot(y[:N], y[N:], step)


create_gif()


torch.cuda.synchronize()
end_time = time.time()
print(f"Elapsed: {end_time - start_time:.2f} s")