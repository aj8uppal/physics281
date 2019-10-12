import numpy as np
import sys, os, time
from numpy.random import RandomState
import matplotlib.pyplot as pl

TAU = float(input("Total optical depth (1.1): ") or 1.1)
NUM_OF_SLICES = int(input("Number of slices (50): ") or 50)
EFFECTIVE_TEMP = int(input("Effective temperature (253): ") or 253)
NUM_OF_SIMUL = int(input("Number of photons (100,000): ") or 100000)
DTAU = TAU/NUM_OF_SLICES

r = RandomState()
layers = np.linspace(0, TAU, NUM_OF_SLICES)
layer_status = [0]*(NUM_OF_SLICES) #The optical depth is defined to be 0 at the top of the atmosphere and increasing with depth into the atmosphere.
photons = {"escaped": 0, "grounded": 0}

def getGroundTemperature():
    ratio = photons["grounded"]/photons["escaped"]
    return (EFFECTIVE_TEMP**4*(1+ratio))**0.25

def getTemperatureOfLayer(i):
    return ((layer_status[i]*EFFECTIVE_TEMP**4)/(2*DTAU*photons["escaped"]))**0.25

def loop(layer_status):
    layer = NUM_OF_SLICES-1
    done = False
    direction = -1
    while (0 <= layer < NUM_OF_SLICES):
        if r.rand() < DTAU:
            layer_status[layer]+=1
            direction = [-1, 1][round(r.rand())]
        layer += direction
    photons["escaped" if layer == -1 else "grounded"]+=1
    return layer_status

print("\nPerforming {} simulations of random walks for photons through {} slices of total optical depth {}...\n".format(NUM_OF_SIMUL, NUM_OF_SLICES, TAU))
time.sleep(0.5)

startTime = time.time()
rows, columns = os.popen('stty size', 'r').read().split()
for i in range(NUM_OF_SIMUL):
    layer_status = loop(layer_status)
    elapsed = round(time.time() - startTime)
    total = round((time.time() - startTime) * (NUM_OF_SIMUL/(i+1)))
    remaining = total - elapsed
    prefix = "\r{}/{} ({}%) [Elapsed time: {} Remaining time: {}] ".format(i+1, NUM_OF_SIMUL, round(i/NUM_OF_SIMUL*100), elapsed, remaining)
    progressbar_len = int(columns) - (len(prefix) + 2)
    progressbar = "".join(["#" if j <= i*progressbar_len//NUM_OF_SIMUL else "." for j in range(progressbar_len)])
    sys.stdout.write(prefix+"["+progressbar+"]");
    sys.stdout.flush()


# print(layer_status)

ratio = photons["grounded"]/photons["escaped"]
groundTemp = getGroundTemperature()
print("\nOut of {} photons, {} were absorbed by the ground and {} were reflected into space. This results in a ratio P_g/P_e of {}.".format(NUM_OF_SIMUL, photons["grounded"], photons["escaped"], round(ratio, 3)))
print("We can calculate the ground temperature by taking the 4th root of T_eff**4*(1+P_g/P_e) = ({0}**4*(1+{1}))**0.25 = {2:.2f}".format(EFFECTIVE_TEMP, round(ratio, 3), groundTemp))
print("For a tau of {0}, and an effective temperature of {1}, the simulated ground temperature is {2:.2f}.".format(TAU, EFFECTIVE_TEMP, getGroundTemperature()))

x_vals = [getTemperatureOfLayer(i) for i in range(0, NUM_OF_SLICES)]
y_vals = range(len(layer_status))

z = np.polyfit(x_vals, y_vals, 1)
f = np.poly1d(z)

x_new = np.linspace(x_vals[0], x_vals[-1], NUM_OF_SLICES)
y_new = f(x_new)

slope, intercept = z

print("y = {}x - {}".format(slope, intercept*-1))

pl.plot(x_new, y_new, color="orange")

print("The top of the atmosphere (layer 1) has a temperature {}.".format(x_vals[0]))

pl.xlim([round(min(x_vals)*0.9524), 300])
pl.ylim([-NUM_OF_SLICES/20, round(NUM_OF_SLICES*1.1)])
pl.xlabel("Temperature (K)")
pl.ylabel(r"Optical Depth ($d\tau$) ")

pl.title("Temperature vs Optical Depth")

plot = pl.scatter(x_vals, y_vals, s=7.5, c='#9FBCFF')
pl.scatter(groundTemp, NUM_OF_SLICES, s=30, c='#002E98')
pl.show()
