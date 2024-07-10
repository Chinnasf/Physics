import numpy as np

# C_{19}HF_{39}O_{5}S >>> rep 1
# C_{38}H_{2}F_{78}O_{10}S_{2} >>> rep 2

elements  = ["C", "H", "F", "O", "S"]
atoms_nmr = np.array([19, 1, 39, 5, 1])

rep  = 2
text = f"/polymer{2**rep}///AQQS`0/"

for i, element in enumerate(elements):
	for j in range( 1, (2**(rep-1))*(atoms_nmr[i])+1 ):
		print(f'alter {text}{elements[i]}{j}, name="{elements[i]}{(atoms_nmr[i]*2**(rep-1))+j}"')

"""
create polymer4, polymer2

# Distance from center of mass to C32 (polymer2)
# 18.156854315862592

translate [18.15*2, 0, 0], polymer4

p2C32 = np.array(cmd.get_atom_coords("/polymer2///AQQS`0/C32"))
p4C17 = np.array(cmd.get_atom_coords("/polymer4///AQQS`0/C17"))

print(np.linalg.norm(p2C32-p4C17))


[18.15*2, -2, -2]: 4.118356345154373
[18.15*2, 0, 0]  : 3.0817572229623584
[17.2*2, 0, 0]   : 1.7525490384384725
[17.05*2, 0, 0]  : 1.6482789868746706

[18*2, -1, -2]: 3.8525817026913827
[18*2, -1, 0] : 2.5742002715393393
[17.6*2, -0.7, 1.5]: 1.6768528786692907



remove /monomer2///AQQS`0/F41
remove /polymer4///AQQS`0/F30


 You clicked /polymer2///AQQS`0/F41
 Selector: selection "sele" defined with 1 atoms.
 You clicked /polymer4///AQQS`0/F31
 Selector: selection "sele" defined with 2 atoms.

"""