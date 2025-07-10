import torch
import os
import numpy as np
import h5py

# Device setup
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')



######### Parameters
L = 100       # Domain of the solution 0 <= x <= L  (in Debye lengths)
N = 200000    # Number of electrons
J = 2**10     # Number of grid-points
vb = 150      # Beam velocity
dt = 0.1      # time step  (in inverse plasma frequencies)
t_max = 60    # such that 0 <= t <= t_max

n0 = N / L    # ion number density
dx = L / J   
timesteps = int(t_max / dt)

# Check input parameters make sence:
if (N < 1) | (J < 2) | (L <= 0.) | (vb <= 0.) | (dt <= 0.) | (t_max <= 0.) | ((int(t_max / dt) / 10) < 1):
    print("Error - invalid input parameters")
    
# To save images
sim_params = f"L{L}_N{N}_J{J}_vb{vb}_dt{dt}_tmax{t_max}_"
os.makedirs("Results", exist_ok=True)

# Preallocated Buffers
ne  = torch.zeros(J, device=device)                         # electron number density
rho = torch.zeros(J, device=device)                         # n(x)/n_0 - 1
rho_hat = torch.zeros(J, dtype=torch.cfloat, device=device) # FFT(rho)
k   = 2*np.pi*torch.fft.fftfreq(J, d=L/J).to(device)        # 2π/L 
phi_hat = torch.zeros(J, dtype=torch.cfloat, device=device) # FFT(electric potential)
phi = torch.zeros(J, device=device)                         # electric potential
E_grid  = torch.zeros(J, device=device)                     # electric field 


# Sample Initial Velocities | Method: rejection sampling
def sample_velocity(v_b, tails_factor, batch_size=100000):
    """
    batch_size is here to avoid generating one sample per loop iteration.
    """
    f_max = 0.5*(1.0 + np.exp(-2*v_b**2))
    vmin, vmax = -tails_factor*v_b, tails_factor*v_b
    velocities = torch.empty(N, device=device)
    filled = 0
    while filled < N: 
        v_ = (vmax - vmin) * torch.rand(batch_size, device=device) + vmin
        beam_shape = 0.5 * (torch.exp(-0.5 * (v_ - v_b)**2) + torch.exp(-0.5 * (v_ + v_b)**2))
        gamma = torch.rand(batch_size, device=device) * f_max
        accepted = v_[gamma <= beam_shape]
        num_to_fill = min(accepted.shape[0], N - filled)
        velocities[filled:filled + num_to_fill] = accepted[:num_to_fill]
        filled += num_to_fill
    return velocities


# Scriptable derivative function
@torch.jit.script
def compute_derivative_jit(y_local: torch.Tensor, N: int, J: int, dx: float, n0: float,
                           ne: torch.Tensor, rho: torch.Tensor,
                           rho_hat: torch.Tensor, phi_hat: torch.Tensor,
                           phi: torch.Tensor, E_grid: torch.Tensor,
                           k: torch.Tensor) -> torch.Tensor:
    r = y_local[:N]
    v = y_local[N:]

    ne.zero_()
    j = torch.floor(r / dx).to(torch.int64)
    f = (r / dx) - j.to(torch.float32)
    ne.scatter_add_(0, (j + 1) % J, f / dx)
    ne.scatter_add_(0, j % J, (1 - f) / dx)

    rho.copy_(ne / n0 - 1)
    rho_hat.copy_(torch.fft.fft(rho))

    phi_hat.zero_()
    phi_hat[1:] = -rho_hat[1:] / (k[1:]**2)
    phi_hat[0] = torch.tensor(0.0 + 0.0j, dtype=torch.cfloat, device=phi_hat.device)
    phi.copy_(torch.fft.ifft(phi_hat).real)
    phi_left = torch.roll(phi, shifts=1, dims=0)
    phi_right = torch.roll(phi, shifts=-1, dims=0)
    E_grid.copy_((phi_left - phi_right) / (2 * dx))

    E_particles = (1 - f) * E_grid[j % J] + f * E_grid[(j + 1) % J]

    return torch.cat([v, -E_particles])



# Improved RK4 using JIT derivative
def RungeKutta4(y_):
    k1 = compute_derivative_jit(y_, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
    k2 = compute_derivative_jit(y_ + 0.5*dt*k1, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
    k3 = compute_derivative_jit(y_ + 0.5*dt*k2, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
    k4 = compute_derivative_jit(y_ +     dt*k3, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
    return y_ + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def run_simulation():
    
    # Initial state
    r0 = torch.rand(N, device=device)*L
    v0 = sample_velocity(vb, 4.0)
    y = torch.cat([r0, v0])
    
    with h5py.File(f"Results/data/{sim_params}.h5", 'w') as h5file:
        h5file.create_dataset("initial_positions", data=r0.cpu().numpy(), dtype='f4')
        h5file.create_dataset("initial_velocities", data=v0.cpu().numpy(), dtype='f4')
        h5file.create_dataset("positions", shape=(timesteps, N), dtype='f4')
        h5file.create_dataset("velocities", shape=(timesteps, N), dtype='f4')
        h5file.create_dataset("time", shape=(timesteps,), dtype='f4')
        h5file.create_dataset("kinetic_energy", shape=(timesteps,), dtype='f4')
        h5file.create_dataset("momentum", shape=(timesteps,), dtype='f4')
        

        for step in range(timesteps):
            y = RungeKutta4(y)
            y[:N] = y[:N] % L

            # Save to file
            h5file["positions"][step] = y[:N].cpu().numpy()
            h5file["velocities"][step] = y[N:].cpu().numpy()
            h5file["time"][step] = step * dt

            kinetic_energy = 0.5 * torch.mean(y[N:]**2).item()
            h5file["kinetic_energy"][step] = kinetic_energy
            
            momentum = torch.sum(y[N:]).item()
            h5file["momentum"][step] = momentum

            

run_simulation()

print(f"\nALL GOOD! Data saved ---> Results/data/{sim_params}.h5")
