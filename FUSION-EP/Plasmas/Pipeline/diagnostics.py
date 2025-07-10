import os
import matplotlib.pyplot as plt
import glob
import imageio.v2 as imageio
import h5py

#plt.style.use(['science', 'notebook'])
plt.style.use(['seaborn-v0_8-notebook'])

def get_data(file):
    with h5py.File("Results/data/" + file, "r") as f:
        time = f["time"][:]
        Momentum = f["momentum"][:]
        Kinetic_E = f["kinetic_energy"][:]
        v_t = f["velocities"][:]
        x_t = f["positions"][:]
        x_0 = f["initial_positions"][:]
        v_0 = f["initial_velocities"][:]
    return (time, x_0, v_0, x_t, v_t, Kinetic_E, Momentum)


def get_simulation_params(file):
    L  = int(file.split('_')[0][1:])
    N  = int(file.split('_')[1][1:])
    J  = int(file.split('_')[2][1:])
    vb = int(file.split('_')[3][2:]) 
    dt = float(file.split('_')[4][2:]) 
    tmax = int(file.split('_')[5][4:])
    return (L, N, J, vb, dt, tmax)


def create_sim_images(file): 
    image_location = "Results/Images/" + file.split('.h5')[0] + "/"
    os.makedirs(image_location, exist_ok=True)
    
    time, _, v_0, x_t, v_t, _, _ = get_data(file)
    L, N, _, vb, dt, _ = get_simulation_params(file)
    velocity_tags = abs(v_0 - vb) < abs(v_0 + vb) # False ==> left

    for t in range(len(time)): 
        plt.plot(x_t[t][velocity_tags], v_t[t][velocity_tags],'.', markersize=1, alpha=0.5, color="firebrick")
        plt.plot(x_t[t][~velocity_tags], v_t[t][~velocity_tags],'.', markersize=1, alpha=0.5, color="k")
        plt.xlim(0, L)
        plt.ylim(-10, 10)
        plt.title(f"Phase Space - Time {round(t*dt,2)}\n {N} Particles")
        plt.xlabel("Position")
        plt.ylabel("Velocity")
        plt.grid(True, alpha=0.5)

        plot_filename = image_location + f"frame_{t:05d}.png"
        plt.savefig(plot_filename, dpi=150)
        plt.close()
        

def create_sim_GIF(file, fps=15): 
    gif_location = "Results/GIFs/" 
    image_location = "Results/Images/" + file.split('.h5')[0] + "/"
    os.makedirs(gif_location, exist_ok=True)
    frames = sorted(glob.glob(image_location+"/frame_*.png"))
    images = [imageio.imread(f) for f in frames]
    imageio.mimsave(gif_location+file.split('.h5')[0]+f"_GIF_fps_{fps}.gif", images, fps=fps)



file = input("file: ")
print("\nok! One moment, please.")

create_sim_images(file)
print("\n\nFinished generating images for the GIF")


create_sim_GIF(file)
print("\n\nDONE!")