from ase import io

# Load the entire trajectory file
traj = io.read('thermal_nvt.traj', index=':')

# Extract the last frame
last_frame = traj[-1]

# Save as a PDB file
last_frame.write('last_frame_ionomer.pdb')

# Save as a LAMMPS data file
#last_frame.write('last_frame_ionomer.lmp')
