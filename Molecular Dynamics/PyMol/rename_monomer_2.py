import numpy as np

# C_{19}HF_{39}O_{5}S

element = "C"
total = 19
rep   = 1

text = f"/monomer{rep+1}///AQQS`0/"

for i in range(1,total+1):
	print(f'alter {text}{element}{i}, name="{element}{rep*total + i}"')


# alter /monomer2///AQQS`0/H1, name="H2"
# alter /monomer2///AQQS`0/S1, name="S2"

"""

load AQQS_allatom_optimised_geometry.pdb, monomer1
import numpy as np

create monomer2, monomer1

rotate x, 180, monomer2
translate [17, -2, -2], monomer2

m1C13 = np.array(cmd.get_atom_coords("/monomer1///AQQS`0/C13"))
m2C17 = np.array(cmd.get_atom_coords("/monomer2///AQQS`0/C17"))

print(np.linalg.norm(m1C13-m2C17))


####################################
### IDEAS 

# Eliminate m1F2 and m2F31, the rotate w.r.t. m2C17.
# O podrias rotar 195 en X y volver a optimizar distancia.


# 3.5797632085931945 [18.1, -2, 0]
# 3.4916854267350836 [18.1, -1.7, 0]
# 3.357916637832542  [18, -1.2, 0]
# 2.712061819565996  [18.6, -1.2, -1.2]
# 1.6243842818508853 [17, -2, -2]
"""
