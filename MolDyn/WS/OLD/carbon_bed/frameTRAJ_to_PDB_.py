from ase import io

# Load the entire trajectory file
traj = io.read('zero_C_bed.traj', index=':')

# Extract the last frame
last_frame = traj[0]

# Save as a PDB file
last_frame.write('carbon_bed.pdb')

# Save as a LAMMPS data file
#last_frame.write('last_frame_ionomer.lmp')
