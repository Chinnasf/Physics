import torch
import os
import matplotlib.pyplot as plt
import numpy as np
import time
import glob
from torch.profiler import profile, record_function, ProfilerActivity, tensorboard_trace_handler, schedule

# Device setup
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
torch.cuda.synchronize()
start_time = time.time()

######### Parameters
L = 100
N = 250000
J = 2**10
vb = 6
n0 = N / L
dx = L / J
dt = 0.5
t_max = 60
timesteps = int(t_max / dt)

if (N < 1) | (J < 2) | (L <= 0.) | (vb <= 0.) | (dt <= 0.) | (t_max <= 0.) | ((int(t_max / dt) / 10) < 1):
    print("Error - invalid input parameters")

# Profiler setup
prof = profile(
    activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA],
    schedule=schedule(wait=1, warmup=1, active=5, repeat=1),
    on_trace_ready=tensorboard_trace_handler("./logdir"),
    record_shapes=True,
    profile_memory=True,
    with_stack=True
)
prof.__enter__()

# Sample initial velocities
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
        velocities[filled:filled + num_to_fill] = accepted[:num_to_fill]
        filled += num_to_fill

    return velocities

# Preallocated buffers
ne = torch.zeros(J, device=device)
rho = torch.zeros(J, device=device)
rho_hat = torch.zeros(J, dtype=torch.cfloat, device=device)
k = 2 * np.pi * torch.fft.fftfreq(J, d=L/J).to(device)
phi_hat = torch.zeros(J, dtype=torch.cfloat, device=device)
phi = torch.zeros(J, device=device)
E_grid = torch.zeros(J, device=device)

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
    with record_function("RungeKutta4"):
        k1 = compute_derivative_jit(y_, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
        k2 = compute_derivative_jit(y_ + 0.5 * dt * k1, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
        k3 = compute_derivative_jit(y_ + 0.5 * dt * k2, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
        k4 = compute_derivative_jit(y_ + dt * k3, N, J, dx, n0, ne, rho, rho_hat, phi_hat, phi, E_grid, k)
        return y_ + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

# Initial state
r0 = torch.rand(N, device=device) * L
v0 = sample_velocity(vb, 4)
y = torch.cat([r0, v0])

# MAIN LOOP
for step in range(timesteps):
    with record_function("main_simulation_step"):
        y = RungeKutta4(y)
        y[:N] = y[:N] % L
    prof.step()

prof.__exit__(None, None, None)
torch.cuda.synchronize()
end_time = time.time()
print(f"Elapsed: {end_time - start_time:.2f} s")
print(prof.key_averages().table(sort_by="cuda_time_total", row_limit=25))
