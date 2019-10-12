import numpy as np
import sys, os, time
from numpy.random import RandomState
import matplotlib.pyplot as pl
from scipy.interpolate import make_interp_spline, BSpline

# TAU = float(input("Total optical depth (1.1): ") or 1.1)

H = 8000 #meters
atmosphere_height = 100000 #100km, this is where "space" begins

def getPositionsOfLayers(target, slices):
    #slices = # of slices
    #target = dm, target mass per layer
    layer_nums = range(slices)
    return list(reversed([-H*np.log(1-((target*layer)/(density(0)*H))) for layer in layer_nums]+[atmosphere_height]))

def density(altitude):
    # p(0)*e^(-altitude/H)
    ground_density = 1.225 #kg/m^2
    return ground_density*np.e**(-altitude/H)

def mass(a0, af):
    def massUnder(altitude):
        return density(0)*H*(1-np.e**(-altitude/H))
    return massUnder(af)-massUnder(a0)


M = density(0)*H
TAU = float(input("Tau (1.1): ") or 1.1)
NUM_OF_SLICES = int(input("Number of slices (50): ") or 50)
EFFECTIVE_TEMP = int(input("Effective temperature (253): ") or 253)
NUM_OF_SIMUL = int(input("Number of photons (100,000): ") or 100000)
DM = M/NUM_OF_SLICES

#what fraction of tau produces a mass of DM
#fraction is probability of absorption

positions = getPositionsOfLayers(DM, NUM_OF_SLICES)
# print(positions)



# for i in range(len(positions)):


r = RandomState()

# layers = np.linspace(0, TAU, NUM_OF_SLICES)
layers = positions
layer_status = [0]*(NUM_OF_SLICES) #The optical depth is defined to be 0 at the top of the atmosphere and increasing with depth into the atmosphere.
photons = {"escaped": 0, "grounded": 0}

def getDTau(i):
    return 1.1*((layers[i]-layers[i+1])/100000)

def getProbabilityFromWidth(i):
    return 1-np.e**(-getDTau(i))

def getGroundTemperature():
    ratio = photons["grounded"]/photons["escaped"]
    return (EFFECTIVE_TEMP**4*(1+ratio))**0.25

def getTemperatureOfLayer(i):
    return ((layer_status[i]*EFFECTIVE_TEMP**4)/(2*getDTau(i)*photons["escaped"]))**0.25

def loop(layer_status):
    layer = NUM_OF_SLICES-1
    done = False
    direction = -1
    while (0 <= layer < NUM_OF_SLICES):
        if r.rand() < getProbabilityFromWidth(layer):
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

altitudes = list(reversed(positions[1:]))[:-1]
y_vals = list(reversed([getTemperatureOfLayer(i) for i in range(0, NUM_OF_SLICES)]))[:-1]

z = np.polyfit(altitudes, y_vals, 3)
f = np.poly1d(z)

# slope, intercept = z

x_new = np.linspace(0, 300, len(altitudes))
y_new = f(x_new)


# pl.plot(x_new, y_new)
# pl.xlim([altitudes[0]-1, altitudes[-1] + 1 ])
# plt.show()
# x_vals = np.linspace(0, 100000, 1000)
# y_vals = range(len(layer_status))

# print("The top of the atmosphere (layer 1) has a temperature {}.".format(x_vals[0]))

# pl.xlim([round(min(altitudes)*0.9524), 300])
# pl.ylim([-NUM_OF_SLICES/20, round(NUM_OF_SLICES*1.1)])
pl.xlabel("Altitude (m)")
pl.ylabel(r"Temperature (K)")

pl.title("Altitude (m) vs Temperature (K)")

plot = pl.scatter(altitudes, y_vals, s=7.5, c='#9FBCFF')
pl.scatter(0, groundTemp, s=30, c='#002E98')
# xnew = np.linspace(0, 30000, 300)
# spl_smooth = make_interp_spline(altitudes, y_vals, k=3) #BSpline object
# smooth = spl_smooth(xnew)
#
# pl.plot(xnew, smooth)

pl.show()
