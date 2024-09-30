import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')

from IPython.core.debugger import Pdb #Pdb().set_trace()


import re
import os

filename_ = "__.data"
output_filename_ = "__.mol"

slices_structure = [ , ]
slices_bonds = [ , ]
slices_angles = [ , ]
slices_dihedral = [ , ]

def get_data(filename_):
	with open(filename_, "r") as file:
		file_content = file.readlines()
	data = pd.DataFrame(file_content, columns=['data'])    
	return data


#df = get_data(filename_)
#Pdb().set_trace()

def get_structure(slices):
	data  = get_data(filename_)
	coords = data.loc[slices[0]:slices[-1]].copy().reset_index(drop=True)
	coords["atom id"] = coords["data"].str.split('\s').apply(lambda x: x[0]).astype(int)
	coords["molecule id"] = coords["data"].str.split('\s').apply(lambda x: x[1]).astype(int)
	coords["atom type"] = coords["data"].str.split('\s').apply(lambda x: x[2]).astype(int)
	coords["charge"] = coords["data"].str.split('\s').apply(lambda x: x[3]).astype(float)
	coords["x"] = coords["data"].str.split('\s').apply(lambda x: x[4]).astype(float)
	coords["y"] = coords["data"].str.split('\s').apply(lambda x: x[5]).astype(float)
	coords["z"] = coords["data"].str.split('\s').apply(lambda x: x[6]).astype(float)
	coords = coords.sort_values("atom id").reset_index(drop=True)
	return coords

def get_bonds(slices):
	data  = get_data(filename_)
	bonds = data.loc[slices[0]:slices[-1]].copy().reset_index(drop=True)
	bonds["bond id"] = bonds["data"].str.split('\s').apply(lambda x: x[0]).astype(int)
	bonds["bond type"] = bonds["data"].str.split('\s').apply(lambda x: x[1]).astype(int)
	bonds["atom1"] = bonds["data"].str.split('\s').apply(lambda x: x[2]).astype(int)
	bonds["atom2"] = bonds["data"].str.split('\s').apply(lambda x: x[3]).astype(int)
	return bonds

def get_angles(slices):
	data  = get_data(filename_)
	angles = data[slices[0]:slices[1]].copy().reset_index(drop=True)
	angles["data"] = angles["data"].str.split('\s').apply(lambda x: list(filter(None, x)))
	angles["angle id"] = angles["data"].apply(lambda x: x[0]).astype(int)
	angles["angle type"] = angles["data"].apply(lambda x: x[1]).astype(int)
	angles["atom1"] = angles["data"].apply(lambda x: x[2]).astype(int)
	angles["atom2"] = angles["data"].apply(lambda x: x[3]).astype(int)
	angles["atom3"] = angles["data"].apply(lambda x: x[4]).astype(int)
	return angles


def get_dihedrals(slices):
	data  = get_data(filename_)
	dihedrals = data[slices[0]:slices[1]].copy().reset_index(drop=True)
	dihedrals["data"] = dihedrals["data"].str.split('\s').apply(lambda x: list(filter(None, x)))
	dihedrals["dihedral id"] = dihedrals["data"].apply(lambda x: x[0]).astype(int)
	dihedrals["dihedral type"] = dihedrals["data"].apply(lambda x: x[1]).astype(int)
	dihedrals["atom1"] = dihedrals["data"].apply(lambda x: x[2]).astype(int)
	dihedrals["atom2"] = dihedrals["data"].apply(lambda x: x[3]).astype(int)
	dihedrals["atom3"] = dihedrals["data"].apply(lambda x: x[4]).astype(int)
	dihedrals["atom4"] = dihedrals["data"].apply(lambda x: x[5]).astype(int)
	return dihedrals


structure = get_structure(slices_structure)
bonds = get_bonds(slices_bonds)
angles = get_angles(slices_angles)
dihedrals = get_dihedrals(slices_dihedral)

#Pdb().set_trace()


def write_mol_file():
    f = open(output_filename_, "w")
    f.write("# Created by Karina Chiñas Fuentes | Python | 26.09.24\n\n")
    f.write(f"\n{len(structure)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n{len(dihedrals)} dihedrals")
    f.write("\n\n")
    f.write("Coords\n\n")
    for i in range(len(structure)):
        f.write(f"{structure['atom id'][i]}\t{structure['x'][i]}\t{structure['y'][i]}\t{structure['z'][i]}\n")
    f.write("\nTypes\n\n")
    for i in range(len(structure)):
        f.write(f"{structure['atom id'][i]}\t{structure['atom type'][i]}\n")
    f.write("\nCharges\n\n")
    for i in range(len(structure)):
        f.write(f"{structure['atom id'][i]}\t{structure['charge'][i]}\n")
    f.write("\nBonds\n\n")
    for i in range(len(bonds)):
        f.write(f"{bonds['bond id'][i]}\t{bonds['bond type'][i]}\t"+
                f"{bonds['atom1'][i]}\t{bonds['atom2'][i]}\n")
    f.write("\nAngles\n\n")
    for i in range(len(angles)):
        f.write(f"{angles['angle id'][i]}\t{angles['angle type'][i]}\t"+
                f"{angles['atom1'][i]}\t{angles['atom2'][i]}\t{angles['atom3'][i]}\n")
    f.write("\nDihedrals\n\n")
    for i in range(len(dihedrals)):
        f.write(f"{dihedrals['dihedral id'][i]}\t{dihedrals['dihedral type'][i]}\t{dihedrals['atom1'][i]}\t"+
                f"{dihedrals['atom2'][i]}\t{dihedrals['atom3'][i]}\t{dihedrals['atom4'][i]}\n")

write_mol_file()

print(f"Created <<{output_filename_}>> file successfully!")