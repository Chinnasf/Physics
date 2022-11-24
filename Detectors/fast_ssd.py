import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import scipy.constants as spc
import pylab as pl
import pandas as pd



Rdetect = 60       # mm
pos_detector = 70  # mm

def compute_solid_angle(N, Rsource, sseedd):
	np.random.seed(sseedd)

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

	# Finding the landing position as in the detector
	x = pos_detector*np.tan(θ)*np.cos(φ) + x0
	y = pos_detector*np.tan(θ)*np.sin(φ) + y0


	r = np.sqrt(np.square(x) + np.square(y))
	hit_index = np.where(r<Rdetect)[0]
	Nhit = len(hit_index)

	Ntot = 2*N
	Ω = 4*np.pi*Nhit/Ntot

	return Ω


# compute_solid_angle(N, Rsource, sseedd)
#compute_solid_angle(975000, 40, 23478)	


r = np.linspace(0.01, 40)
for i in r:
	print(f"{i},{compute_solid_angle(int(1e7), i, 23478)}")