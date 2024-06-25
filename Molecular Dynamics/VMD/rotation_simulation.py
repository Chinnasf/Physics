import numpy as np

angle = np.linspace(0,360,50)

print("display resetview")
for i in np.linspace(0,360,50):
	angle = round(i, 2)
	print(f"render snapshot frame_{angle}.rgb")
	print("display resetview")
	print(f"rotate y by {angle}")
	