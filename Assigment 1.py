
# Air density
P = 1e5
R = 287
T = 25 + 273.15

rho = P / (R * T)   # Air density in kg / m3

import numpy as np
import matplotlib.pyplot as plt
import efficiencyCalculations
import trajectoryCalculations
#this is just filled with dummy values for now
#Marcelo please update this to have points along the boundary
boundary_matrix = np.array([
    [0,0],
    [0,1]
])


def cot(x):
    return 1 / np.tan(x)

# Range of t (in radians)
t = np.linspace(np.deg2rad(15.5568), np.deg2rad(176.891), 500)

# Curve 1
C1 = 3.7699111843077517
x1 = (5*(C1 - t) * cot(t)) / (2*np.pi)
y1 = 5*(C1 - t) / (2*np.pi)

# Curve 2
C2 = 4.39822971502571
x2 = (5*(C2 - t) * cot(t)) / (2*np.pi)
y2 = -5*(C2 - t) / (2*np.pi)

# Curve 3
x3 = (5*(np.pi - t) * cot(t)) / (2*np.pi)
y3 = -5*(np.pi - t) / (2*np.pi)

# Curve 4
x4 = (5*(np.pi - t) * cot(t)) / (2*np.pi)
y4 = 5*(np.pi - t) / (2*np.pi)

#Particle Trajectory -------------------------------------------------------------
trajectoryCalculations.update_boundary_matrix(boundary_matrix)
#calling with dummy values for now
test_matrix = trajectoryCalculations.run_particle_simulation(startPos=[0,0],startVelocity=[2,4],startAirFlow=[1,2])

# Define integration range only for Inlet
mask_inlet = (x1 >= -10) & (x1 <= -0.8)
# Define integration range for Scavenger (Outlet)
mask_scavenger = (x1 > -0.8) & (x1 <= 10)
# Define integration range for Outlet Core
mask_core = (x2 > -0.8) & (x2 <= 10)

# Inlet area
x_inlet = x1[mask_inlet]
y_inlet = y1[mask_inlet] - y2[mask_inlet]

# Sort if necessary
if x_inlet[0] > x_inlet[-1]:
    x_inlet = x_inlet[::-1]
    y_inlet = y_inlet[::-1]

A_inlet = np.trapezoid(y_inlet, x_inlet)
print(f"Inlet area = {A_inlet:.2f} m²")

# Select x values and height for Scavenger
x_scavenger = x1[mask_scavenger]
y_scavenger = y1[mask_scavenger] - y4[mask_scavenger]

# Sort if necessary
if x_scavenger[0] > x_scavenger[-1]:
    x_scavenger = x_scavenger[::-1]
    y_scavenger = y_scavenger[::-1]

# Integrate using trapezoid rule
A_scavenger = np.trapezoid(y_scavenger, x_scavenger)
print(f"Scavenger area = {A_scavenger:.2f} m²")

# Select x values and height for Core
x_core = x2[mask_core]
y_core = y3[mask_core] - y2[mask_core]

# Sort if necessary
if x_core[0] > x_core[-1]:
    x_core = x_core[::-1]
    y_core = y_core[::-1]

# Integrate using trapezoid rule
A_core = np.trapezoid(y_core, x_core)
print(f"Core area = {A_core:.2f} m²")


V_inlet = 1.0  # m/s
m_dot_inlet = rho * V_inlet * A_inlet
print(f"Mass flow rate at inlet = {m_dot_inlet:.3f} kg/s")

#post simmulation efficiency calculation
#need to pass array of masses in each outlet here using dummy for now
core = [5,10,1,5,8,10]
scavenge = [6,4,7,10,1,5,7]
ce = efficiencyCalculations.captureEfficiency(core, scavenge)
print(f"Capture Efficency = {ce:.3f}")



# Plot
plt.figure(figsize=(10,6))
plt.plot(x1, y1, 'b', label='Curve I')
plt.plot(x2, y2, 'g', label='Curve II')
plt.plot(x3, y3, 'r', label='Curve III')
plt.plot(x4, y4, 'orange', label='Curve IV')

# Fill only the Inlet region
plt.fill_between(x1[mask_inlet], y1[mask_inlet], y2[mask_inlet], color='skyblue', alpha=0.4, label='Inlet')
# Fill the Scavenger region
plt.fill_between(x1[mask_scavenger], y1[mask_scavenger], y4[mask_scavenger], color='lightgreen', alpha=0.4, label='Scavenger')
# Fill the Outlet Core region
plt.fill_between(x2[mask_core], y2[mask_core], y3[mask_core], color='salmon', alpha=0.4, label='Core')

plt.axvline(0, color='k', linestyle='--', linewidth=1)
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Inertia particle separator")
plt.grid(True)

# Set axis limits
plt.xlim(-10, 10)
plt.ylim(-10, 10)

plt.show()







