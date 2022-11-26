"""
Fast Solid State Detector code
By: Karina Chinas Fuentes
Student Nummer 02118434
Nov 26, 2022
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import scipy.constants as spc
import pylab as pl
import seaborn as sns
import pandas as pd
import time


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

sstyle = "paper"
sns.set_context(sstyle)
sns.set_style("whitegrid")
plt.rc('font',family = 'serif')

Rdetect = 60       # mm
pos_detector = 70  # mm

def compute_x_y_pos(N, Rsource, sseedd):
	np.random.seed(sseedd)

	N = int(N)

	# Creating values for the position
	x_ = np.random.uniform(-Rsource, Rsource, N)
	y_ = np.random.uniform(-Rsource, Rsource, N)

	r_ = np.sqrt( np.square(x_) + np.square(y_) )
	index_ = np.where(r_ < Rsource)

	y0 = y_[index_]
	x0 = x_[index_]
	N = len(index_[0])

	# Creating values for the direction
	θ = np.arccos(2*np.random.rand(N) - 1)
	φ = np.pi*np.random.rand(N)

	# Finding the landing position as in the detector pos
	x = pos_detector*np.tan(θ)*np.cos(φ) + x0
	y = pos_detector*np.tan(θ)*np.sin(φ) + y0

	r = np.sqrt(np.square(x) + np.square(y))
	hit_index = np.where(r<Rdetect)[0]
	Nhit = len(hit_index)

	Ntot = 2*N
	Ω = 4*np.pi*Nhit/Ntot

	return x,y



# Solution to problem 1 
def compute_solid_angle(N, Rsource, sseedd):
	x,y = compute_x_y_pos(N, Rsource, sseedd)	
	r = np.sqrt(np.square(x) + np.square(y))
	hit_index = np.where(r<Rdetect)[0]
	Nhit = len(hit_index)
	Ntot = 2*len(x)
	Ω = 4*np.pi*Nhit/Ntot
	return Ω
#print(compute_solid_angle(1e7, 40, 71))



# Solution to problem 2
# int(1e7), 23478 -- plot
def compute_sa_point_src(N, sseedd):
	r = np.linspace(0.01, 40)
	print("r,omega")
	for i in r:
		print(f"{i},{compute_solid_angle(N, i, sseedd)}")
	return 0
#compute_sa_point_src(1e7, 71)	



# Solution to problem 3
def compute_sqr_opening(N, Rsource, sseedd, percentace = 0.75):
	x,y = compute_x_y_pos(N, Rsource, sseedd)	

	L = (percentace*Rdetect/2)*(2**-0.5)

	df = pd.DataFrame(x,columns=["x"])
	df["y"] = y
	df = df[(df.x > -L) & (df.x < L)]
	df = df[(df.y > -L) & (df.y < L)]
	Nsq_hit = len(df)

	Ntot = 2*len(x)
	Ω = 4*np.pi*Nsq_hit/Ntot
	return 0
#compute_sqr_opening(1e7, 40, 71, percentace = 0.75)



# Solution to problem 4
def plot_fluctuations(Rsource):
	seeds = np.random.randint(10000, size=50)
	seeds = np.unique(seeds)

	Ns = np.array([1e2,1e3,1e4,1e5,1e6,1e7]).astype(int)
	df = pd.DataFrame(seeds, columns=["seeds"])
	for n in Ns:
	    df[f"{n}"] = df.seeds.apply(lambda x: compute_solid_angle(n, 40, x))
	df_sts = df[df.columns[1:]].describe()
	data = (df_sts.loc["std"]/df_sts.loc["mean"])*100

	plt.scatter(data.keys(), data.values, c="r", label="RSD per N")
	plt.axhline(0.5,c="k", label="RSD = 0.5")
	plt.xlabel("Initial Number of Particles")
	plt.ylabel("Coefficient of Variation [%]")
	plt.legend()
	plt.show()

	Ns = np.linspace(1e4,1e6,7).astype(int)
	df = pd.DataFrame(seeds, columns=["seeds"])
	for n in Ns:
	    df[f"{n}"] = df.seeds.apply(lambda x: compute_solid_angle(n, 40, x))
	df_sts = df[df.columns[1:]].describe()
	data = (df_sts.loc["std"]/df_sts.loc["mean"])*100

	plt.scatter(data.keys(), data.values, c="r", label="RSD per N")
	plt.axhline(0.5,c="k", label="RSD = 0.5")
	plt.xlabel("Initial Number of Particles")
	plt.ylabel("Coefficient of Variation [%]")
	plt.legend()
	plt.show()

	return 0
#plot_fluctuations(40)