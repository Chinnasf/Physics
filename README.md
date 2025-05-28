# Physics

<code style="color : red"> Warning: it may take long before all images display.</code>

This repository contains various codes and documents about pure and applied physics developed during my undergraduate and master studies, as well as my PhD research on PEMFC. I will focus on this repository, mainly applying GPU processing to Physics-related insights. For instance, here is an example of a simulation of gas inside a box. The following image shows the code for 500 particles processed with CPU.  

<div align="center">
  <img src="https://github.com/Chinnasf/Physics/blob/master/GPU%20Code/GIFs/withoutGPU_500_smaller__particles.mp4.gif" width="600">
</div>

The following image was created using the same code, but instead of processing with `numpy`, I used `PyTorch` to use the GPU; now, for 5000 particles.

<div align="center">
  <img src="https://github.com/Chinnasf/Physics/blob/master/GPU%20Code/GIFs/withGPU_5000_smaller__particles______.gif" width="600">
</div>

You can also find the code to simulate the damped harmonic/anharmonic oscillator. The code is based on the [Stochastic Processes: Data Analysis and Computer Simulation](https://learning.edx.org/course/course-v1:KyotoUx+009x+1T2017/home) course. 

<div align="center">
  <img src="https://github.com/Chinnasf/Physics/blob/master/General%20Physics/Gifs/AHO_zeta0_005.gif" width="400" align="center">
</div>

### MY WORK AS A MSCA PhD RESEARCHER (DC15) | PROTON EXCHANGE MEMBRANE FUEL CELLS.

The DC15 position is responsible for optimizing the catalyst layer (CL) of a Proton Exchange Membrane Fuel
Cell (PEMFC). There are various types of CLs, with the most common ones relying on platinum as a primary
catalytic component. This type of CL can account for up to 40% of the total cost of a PEMFC [[5]](https://www.sciencedirect.com/science/article/pii/S2451910321001617).
CLs also depend on a polymer to optimize proton transport; hence, the polymer directly impacts both the
performance and longevity of the fuel cell. The most commonly chosen polymers are PFSA-based ionomers,
as they offer commercial advantages due to their lower cost—compared to alternative materials—and their
large-scale production availability [[6]](https://www.sciencedirect.com/science/article/pii/S2214993723001628).

However, increasing research efforts are focused on investigating alternative materials for both the catalyst and
ionomers, aiming to reduce costs and minimize the environmental impact of these devices [[5]](https://www.sciencedirect.com/science/article/pii/S2451910321001617),[[7]](https://www.sciencedirect.com/science/article/pii/S2041652024008435). Although
PFSA-based ionomers, such as Nafion [[8]](https://www.sciencedirect.com/science/article/pii/S0360128510000511), demonstrate excellent engineering performance, they pose significant
risks to both human health and the environment. Due to these concerns, the European Commission has
implemented restrictions on the use of a specific subgroup of PFAS chemicals [[9]](https://ec.europa.eu/commission/presscorner/detail/en/ip_24_4763), [[10]](https://echa.europa.eu/hot-topics/perfluoroalkyl-chemicals-pfas).

Therefore, my position focuses on investigating and optimizing various CL-related materials under different
operating conditions while considering their performance, cost, and environmental impact. To achieve this, the
DC15 is expected to generate its own dataset from molecular dynamics simulations—designed to be as realistic
as possible—which will then be used to train a Physics-Informed Machine Learning algorithm [[3]](https://euraxess.ec.europa.eu/jobs/111519).

The following .gif file results from a simulation I created using LAMMPS. The parametrization of the components were done through a modified version of the  [DREDIGING Force Field](https://pubs.acs.org/doi/10.1021/j100389a010). 

<div align="center">
  <img src="https://github.com/Chinnasf/Physics/blob/master/MolDyn/Simulations/CL_4IL13p5_1p2ns_zoom_in.gif" width="400" align="center">
</div>
---

### Content: OLD Description

- **MolecularModeling**
  - `Ag_simulation.c`: C-code on molecular dynamics code simulating 64 particles of Ag inside a box. The code considers newtonian dynamics alongsine classical statistical mechanics. The code computes the evolution of kinetic energy and potential energy. Libraries used:
    - `<stdio.h>`
    - `<math.h>`
  - `LAMMPS_reader.py`: Python-code used to plot the results of simulations created through [LAMMPS](https://lammps.sandia.gov/). Libraries used:
    - `numpy`
    - `pandas`
    - `matplotlib`
    - `os`
    
  Everything else are visualizations and reports on processes regarding these topics and softwares.
  
- **QuantumMechanics**
  - `EnergyQubits_for_QuantumComputer.py`: Python-code developed during research stay on Quantum Computing. The code computes the energy per state given the number of qubits. The computer is based on a diamond with C13 particles (used as qubits, since they have 1/2 spin). Libraries used:
    - `numpy`
    - `matplotlib`

  - `QuantumTeleportation_4qubits.f95`: Fortran95-code developed during research stay on Quantum Computing. The code simulates a quantum teleportation phenomena for 4 qubits (based on the same system of diamond with C13 particles). 
  
  Everything else are reports on processes regarding quantum mechanics alongside statistical mechanics.
  
- **XPSanalysis**<br>
  - Depth Profile Report: report on a depth profile for a silicon dioxide wafer conducted with an [X-Ray photoelectron spectroscopy](https://www.sciencedirect.com/topics/chemistry/x-ray-photoelectron-spectroscopy) procedure at [ITESO's nanotechnology laboratoies](https://iteso.mx/web/general/detalle?group_id=6542872).
  
  - [Elements Detected in thin film by XPS.pdf](https://github.com/Chinnasf/Physics/blob/master/XPSanalysis/Elements%20Detected%20in%20thin%20film%20by%20XPS.pdf) is an example of detected elements in a thin film.
  
  All visualizations (in pdf format) were created by me using the XPS's data and a Python code similar to `LAMMPS_reader.py`. 
  


---

### Contributing

When contributing to this repository, please first discuss the change you wish to make via email 
(or any other method) with me before making a change.

---

### License

This project is licensed under the [GNU Affero General Public License v3.0](https://www.gnu.org/licenses/agpl-3.0.en.html) -- 
see the [LICENSE](https://github.com/Chinnasf/Physics/blob/master/LICENSE) file for details

