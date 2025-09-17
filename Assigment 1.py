# Air density
P = 1e5
R = 287
T = 25 + 273.15

rho = P / (R * T)   # Air density in kg / m3

import numpy as np
import matplotlib.pyplot as plt

def cot(x):
    return 1 / np.tan(x)

def air_velovity(x,y):
    # Convert to polar coordinates
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)

    # Compute Vr and Vt
    Vr = -np.sin(theta)
    Vt = (1 / (2 * np.pi * r)) + np.cos(theta)

    # Convert to Cartesian velocities
    Vx = Vr * np.cos(theta) - Vt * np.sin(theta)
    Vy = Vr * np.sin(theta) + Vt* np.cos(theta)

    return (Vx, Vy)

# Range of t (in radians)
t = np.linspace(np.deg2rad(15.5568), np.deg2rad(176.891), 5000)

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

# Target x-value
target_x = -10

# Find index where x is closest to target
idx1 = np.argmin(np.abs(x1 - target_x))
idx2 = np.argmin(np.abs(x2 - target_x))

# Get the corresponding points
x1_near, y1_near = x1[idx1], y1[idx1]
x2_near, y2_near = x2[idx2], y2[idx2]

# Air velocity
V_inlet = 1.0  # m/s

# Section length
L_inlet = y1_near - y2_near

# Mass flow rate
m_dot_inlet = rho * V_inlet * L_inlet
print(f"Mass flow rate at inlet = {m_dot_inlet:.4f} kg/s")



# Target x-value for curve 1
target_x = 10
idx = np.argmin(np.abs(x1 - target_x))
x1_scavenger = x1[idx]
y1_scavenger = y1[idx]


# Convert to polar coordinates
r1_scavenger = np.sqrt(x1_scavenger**2 + y1_scavenger**2)
theta1_scavenger = np.arctan2(y1_scavenger, x1_scavenger)
"""
# Compute Vr and Vt
Vr1_scavenger = -np.sin(theta1_scavenger)
Vt1_scavenger = (1/(2*np.pi*r1_scavenger)) + np.cos(theta1_scavenger)

# Convert to Cartesian velocities
Vx1_scavenger = Vr1_scavenger*np.cos(theta1_scavenger) - Vt1_scavenger*np.sin(theta1_scavenger)
Vy1_scavenger = Vr1_scavenger*np.sin(theta1_scavenger) + Vt1_scavenger*np.cos(theta1_scavenger)
"""
Vx1_scavenger, Vy1_scavenger = air_velovity(x1_scavenger,y1_scavenger)


# Target x-value for curve 4
target_x = 10
idx4 = np.argmin(np.abs(x4 - target_x))
x4_scavenger = x4[idx4]
y4_scavenger = y4[idx4]

# Convert to polar coordinates
r4_scavenger = np.sqrt(x4_scavenger**2 + y4_scavenger**2)
theta4_scavenger = np.arctan2(y4_scavenger, x4_scavenger)

# Compute Vr and Vt
Vr4_scavenger = -np.sin(theta4_scavenger)
Vt4_scavenger = (1/(2*np.pi*r4_scavenger)) + np.cos(theta4_scavenger)

# Convert to Cartesian velocities
Vx4_scavenger = Vr4_scavenger*np.cos(theta4_scavenger) - Vt4_scavenger*np.sin(theta4_scavenger)
Vy4_scavenger = Vr4_scavenger*np.sin(theta4_scavenger) + Vt4_scavenger*np.cos(theta4_scavenger)

# Average Y at outlet (Scavenger)
y_avg_scavenger = (y1_scavenger + y4_scavenger) / 2


# Average Vx at outlet (Scavenger)
Vx_avg_scavenger = (Vx1_scavenger + Vx4_scavenger) / 2
Vy_avg_scavenger = (Vy1_scavenger + Vy4_scavenger) / 2

# Outlet section length
L_outlet = y1_scavenger - y4_scavenger

# Mass flow rate at Scavenger outlet
m_dot_scavenger = rho * Vx_avg_scavenger * L_outlet

print(f"\nVelocity at Outlet (Scavenger) At x = {10} m, and y = {y_avg_scavenger:.3f} m: Vx = {Vx_avg_scavenger:.3f} m/s, Vy = {Vy_avg_scavenger:.3f} m/s")

print(f"Mass flow rate at Scavenger outlet = {m_dot_scavenger:.4f} kg/s")



# Target x-value for curve 2
target_x = 10
idx2 = np.argmin(np.abs(x2 - target_x))
x2_core = x2[idx2]
y2_core = y2[idx2]

# Convert to polar coordinates
r2_core = np.sqrt(x2_core**2 + y2_core**2)
theta2_core = np.arctan2(y2_core, x2_core)

# Compute Vr and Vt
Vr2_core = -np.sin(theta2_core)
Vt2_core = (1/(2*np.pi*r2_core)) + np.cos(theta2_core)

# Convert to Cartesian velocities
Vx2_core = Vr2_core*np.cos(theta2_core) - Vt2_core*np.sin(theta2_core)
Vy2_core = Vr2_core*np.sin(theta2_core) + Vt2_core*np.cos(theta2_core)


# Target x-value for curve 3
target_x = 10
idx3 = np.argmin(np.abs(x3 - target_x))
x3_core = x3[idx3]
y3_core = y3[idx3]

# Convert to polar coordinates
r3_core = np.sqrt(x3_core**2 + y3_core**2)
theta3_core = np.arctan2(y3_core, x3_core)

# Compute Vr and Vt
Vr3_core = -np.sin(theta3_core)
Vt3_core = (1/(2*np.pi*r3_core)) + np.cos(theta3_core)

# Convert to Cartesian velocities
Vx3_core = Vr3_core*np.cos(theta3_core) - Vt3_core*np.sin(theta3_core)
Vy3_core = Vr3_core*np.sin(theta3_core) + Vt3_core*np.cos(theta3_core)


# Average Y at outlet (Core)
y_avg_core = (y2_core + y3_core) / 2

# Average Vx at outlet (Core)
Vx_avg_core = (Vx2_core + Vx3_core) / 2
Vy_avg_core = (Vy2_core + Vy3_core) / 2

# Outlet section length (Core)
L_outlet_core = y3_core - y2_core


# Mass flow rate at Core outlet
m_dot_core = rho * Vx_avg_core * L_outlet_core

print(f"\nVelocity at Outlet (Core) At x = {10} m, and y = {y_avg_core:.3f} m: Vx = {Vx_avg_core:.3f} m/s, Vy = {Vy_avg_core:.3f} m/s")

print(f"Mass flow rate at Core outlet = {m_dot_core:.4f} kg/s")
############################################################################################################################################





