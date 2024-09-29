from ase import io

# Load the entire trajectory file
traj = io.read('thermal_nvt___2nd.traj', index=':')

# Extract the last frame
last_frame = traj[169]

# Save as a PDB file
last_frame.write('frame_169_ionomer_nvt.pdb')

# Save as a LAMMPS data file
#last_frame.write('last_frame_ionomer.lmp')
