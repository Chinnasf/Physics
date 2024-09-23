load  ATB_optimized_geometry_3reps.pdb, strc1
import numpy as np


load axes.pml

create strc2, strc1

remove /strc1///Q6MT`0/C50
remove /strc1///Q6MT`0/F43
remove /strc1///Q6MT`0/F44 
remove /strc1///Q6MT`0/F45

remove /strc2///Q6MT`0/C1
remove /strc2///Q6MT`0/FC7
remove /strc2///Q6MT`0/FC8 
remove /strc2///Q6MT`0/FC9

rotate x, -90, strc2
translate [-21.2,-3,16], strc2


STR1C = np.array(cmd.get_atom_coords("/strc1///Q6MT`0/C49"))
STR2C = np.array(cmd.get_atom_coords("/strc2///Q6MT`0/C2"))

# General Reference
C1_MainChain = np.array(cmd.get_atom_coords("/strc1///Q6MT`0/C7"))
C2_MainChain = np.array(cmd.get_atom_coords("/strc1///Q6MT`0/C8"))

print(np.linalg.norm(STR1C-STR2C))
print(" ")
print(np.linalg.norm(C1_MainChain-C2_MainChain))

alter /strc2///Q6MT`0/C2, name="C2unif"
alter /strc1///Q6MT`0/C49, name="C49unif"

create rep1, strc1 or strc2


bond /rep1///Q6MT`0/C2unif, /rep1///Q6MT`0/C49unif

delete strc1
delete strc2

create strc2, rep1

remove /rep1///Q6MT`0/FC7
remove /rep1///Q6MT`0/FC8
remove /rep1///Q6MT`0/FC9
remove /rep1///Q6MT`0/C1

remove /strc2///Q6MT`0/F43
remove /strc2///Q6MT`0/F44
remove /strc2///Q6MT`0/F45
remove /strc2///Q6MT`0/C50


translate [40,0,0], strc2
rotate x, 180, strc2
translate [2.48,15.5,3.35], strc2


zoom animate=1
set mouse_selection_mode = 0


STR1C = np.array(cmd.get_atom_coords("/rep1///Q6MT`0/C2"))
STR2C = np.array(cmd.get_atom_coords("/strc2///Q6MT`0/C49"))

# General Reference
C1_MainChain = np.array(cmd.get_atom_coords("/rep1///Q6MT`0/C7"))
C2_MainChain = np.array(cmd.get_atom_coords("/rep1///Q6MT`0/C8"))

print(np.linalg.norm(STR1C-STR2C))
print(" ")
print(np.linalg.norm(C1_MainChain-C2_MainChain))

alter /rep1///Q6MT`0/C2, name="C2sr"
alter /strc2///Q6MT`0/C49, name="C49sr"

create rep2, rep1 or strc2

bond /rep2///Q6MT`0/C49sr, /rep2///Q6MT`0/C2sr

delete rep1
delete strc2

save polymer_12_reps.pdb


