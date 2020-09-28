import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import rc
import pandas as pd 
import os

"""
DATE: NOV 2019
AUTHOR: Karina Chinas

This code was created to read the results of molecular dynamics simulations
done in LAMMPS. 
"""


sstyle = "seaborn-poster"
plt.style.use(sstyle)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# Open the txt file but it does not recieve it as it is 
# originally written; instead, it is all as a whole 'sentence'. 


choose = 3
#int(input("\n[1]Grafica de Variables Termodinamicas\n[2]Grafica Temp-Densidad y Temp-Energia\n[3]RDF\n\nElige una Opcion: "))
svf = input("\n Do you want to save the figure?[Y/N]: ")


if choose == 1:
	file_name = "R1CSV"

	df = pd.read_csv(file_name + ".csv", sep = ",", header = None, skipinitialspace = True)
	df.columns = ["Step","TempC","E_pair","E_mol","Etot","Press","Volumen"]


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------

	plt.figure(figsize = (20,11))

	plt.plot(df["Step"], df["TempC"], color = "k", linewidth = 0.7)

	#plt.title("Sistema de Agua", fontsize = 25)

	plt.ylabel("Temperatura Teorica (Kelvin)", fontsize = 25)
	plt.xlabel("Tiempo (fs)", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig(file_name + "Temp.png", format = "png")

	plt.show()


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------

	plt.figure(figsize = (20,11))

	plt.plot(df["Step"], df["E_pair"], color = "k", linewidth = 0.7)

	#plt.title("Sistema de Agua", fontsize = 25)

	plt.ylabel("Energia Potencial (Kcal/mol)", fontsize = 25)
	plt.xlabel("Tiempo (fs)", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig(file_name + "Epot.png", format = "png")

	plt.show()


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------

	plt.figure(figsize = (20,11))

	plt.plot(df["Step"], df["Etot"], color = "k", linewidth = 0.7)

	#plt.title("Sistema de Agua", fontsize = 25)

	plt.ylabel("Energia Total (Kcal/mol)", fontsize = 25)
	plt.xlabel("Tiempo (fs)", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig(file_name + "Etot.png", format = "png")

	plt.show()


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------

	plt.figure(figsize = (20,11))

	plt.plot(df["Step"], df["Press"], color = "k", linewidth = 0.7)

	#plt.title("Sistema de Agua", fontsize = 25)

	plt.ylabel("Presion (ATM)", fontsize = 25)
	plt.xlabel("Tiempo (fs)", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig(file_name + "Press.png", format = "png")

	plt.show()


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------

	plt.figure(figsize = (20,11))

	plt.plot(df["Step"], df["Volumen"], color = "k", linewidth = 0.7)

	#plt.title("Sistema de Agua", fontsize = 25)

	plt.ylabel("Volumen (\\AA$^3$)", fontsize = 25)
	plt.xlabel("Tiempo (fs)", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig(file_name + "Vol.png", format = "png")

	plt.show()


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------


elif choose == 2:

	T = [300 + 10*i for i in range(9)]

	Files = [f"RCSV{t}" for t in T]

	DF = [pd.read_csv(Files[i] + ".csv", sep = ",", header = None) for i in range(len(T))]

	for i in range(len(DF)):
		DF[i].columns = ["Step","Temp","E_pair","E_mol","TotEng","Press","Volume"]

	df_E = [DF[i]["TotEng"] for i in range(len(DF))]
	Eprom = []
	for i in range(len(df_E)):
		E = 0
		for j in range(len(DF[0]["TotEng"])):
			E = E + DF[i]["TotEng"][j]
		Eprom.append(E/len(DF[0]["TotEng"]))

	t = np.linspace(min(T),max(T), len(T))
	t = 0.8971359827*t - 31284.941052634

	plt.figure(figsize = (17,10))

	plt.plot(T,t, linewidth = 0.9, label = "Ajuste de Datos", color = 'k')
	plt.scatter(T,Eprom, label = "Datos de LAMMPS",color = 'k')

	plt.text(300,-30950, "$f$(T) = 0.8971359827$\\cdot$T - 31284.941052634", fontsize = 25)

	#plt.title("Sistema de Agua \n Variacion en Temperatura", fontsize = 35)
	plt.ylabel("Energia Total (Kcal/mol)", fontsize = 25)
	plt.xlabel("Temperatura (Kelvin)", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig("Etot_T.png", format = "png")

	plt.show()

	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------

	df_V = [DF[i]["Volume"] for i in range(len(DF))]
	Vprom = []
	for i in range(len(df_V)):
		V = 0
		for j in range(len(DF[0]["Volume"])):
			V = V + DF[i]["Volume"][j]
		Vprom.append(V/len(DF[0]["Volume"]))

	Dnsty = np.zeros(len(Vprom))
	for i in range(len(Vprom)):
		Dnsty[i] = (5**3)*(18.01528/((Vprom[i])*(6.022e23)*(1e-24)))

	t = np.linspace(min(T),max(T), len(T))
	t = 0.000138931833*t + 0.9609146011111

	plt.figure(figsize = (17,10))

	plt.plot(T,t, linewidth = 0.9, label = "Ajuste de Datos", color = 'k')
	plt.scatter(T,Dnsty, label = "Datos de LAMMPS",color = 'k')

	plt.text(300,1.01305, "$f$(T) = 0.000138931833$\\cdot$T + 0.9609146011111", fontsize = 25)

	#plt.title("Sistema de Agua \n Variacion en Temperatura", fontsize = 35)
	plt.ylabel("Densidad (g/cm$^3$)", fontsize = 25)
	plt.xlabel("Temperatura (Kelvin)", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig("Dens_T.png", format = "png")

	plt.show()


elif choose == 3:

	svf = input("\nQuieres guardar las imagenes?[Y/N]: ")

	#T = [240 + 20*i for i in range(7)]
	T = [3000]
	Files = [f"Ni-rdf{t}" for t in T]
	DF = [pd.read_csv(Files[i] + ".dat", sep = " ", header = None) for i in range(len(T))]
	for i in range(len(DF)):
		DF[i].columns = ["temp","r","DRFoo","NCoo","DRFhh","NChh","DRFoh","NCoh"]

	r = DF[0]["r"]
	Goo = [DF[i]["DRFoo"] for i in range(len(T))]


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------


	plt.figure(figsize = (17,10))

	for i in range(len(T)):
		plt.plot(r,Goo[i], label = f"Temperatura = {T[i]} $\\circ$ C")

	#plt.title("Sistema de Agua\n $g_{oo}(r)$", fontsize = 25)
	plt.ylabel("Funcion Radial", fontsize = 25) 
	plt.xlabel("$r$", fontsize = 25) 

	if (svf == "Y" or svf == "y" or svf == ""):
		plt.savefig("NI.png", format = "png")

	plt.legend()	
	plt.show()


	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------



	# --------------------------------------------------------------------------------------
	# --------------------------------------------------------------------------------------




else:
	pass