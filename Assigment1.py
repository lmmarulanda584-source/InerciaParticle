# Air density
import trajectoryCalculations

P = 1e5
R = 287
T = 25 + 273.15

rho = P / (R * T)   # Air density in kg / m3

import numpy as np
import matplotlib.pyplot as plt

def cot(x):
    return 1 / np.tan(x)


# Range of t (in radians)
t = np.linspace(np.deg2rad(15.5568), np.deg2rad(176.891), 5000)

# Curve 1

C1 = 3.7699111843077517
x1 = (5*(C1 - t) * cot(t)) / (2*np.pi)
y1 = 5*(C1 - t) / (2*np.pi)
curve1_points = np.column_stack((x1, y1))

# Curve 2
C2 = 4.39822971502571
x2 = (5*(C2 - t) * cot(t)) / (2*np.pi)
y2 = -5*(C2 - t) / (2*np.pi)
curve2_points = np.column_stack((x2, y2))

# Curve 3
x3 = (5*(np.pi - t) * cot(t)) / (2*np.pi)
y3 = -5*(np.pi - t) / (2*np.pi)
curve3_points = np.column_stack((x3, y3))


# Curve 4
x4 = (5*(np.pi - t) * cot(t)) / (2*np.pi)
y4 = 5*(np.pi - t) / (2*np.pi)
curve4_points = np.column_stack((x4, y4))

# Combine all curves into one matrix of (x,y)
boundary_matrix = np.vstack((curve1_points, curve2_points, curve3_points, curve4_points))
trajectoryCalculations.update_boundary_matrix(boundary_matrix)

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

print(f"\nVelocity at Inlet  At x = {-10} m, and y = {0:.3f} m: Vx = {V_inlet:.3f} m/s")

print(f"Mass flow rate at inlet = {m_dot_inlet:.4f} kg/s")



# Target x-value for curve 1
target_x = 10
idx = np.argmin(np.abs(x1 - target_x))
x1_scavenger = x1[idx]
y1_scavenger = y1[idx]


# Target x-value for curve 4
target_x = 10
idx4 = np.argmin(np.abs(x4 - target_x))
x4_scavenger = x4[idx4]
y4_scavenger = y4[idx4]

# Average Y at outlet (Scavenger)
y_avg_scavenger = (y1_scavenger + y4_scavenger) / 2

# Outlet section length
L_outlet_scavenger = y1_scavenger - y4_scavenger

# Target x-value for curve 2
target_x = 10
idx2 = np.argmin(np.abs(x2 - target_x))
x2_core = x2[idx2]
y2_core = y2[idx2]



# Target x-value for curve 3
target_x = 10
idx3 = np.argmin(np.abs(x3 - target_x))
x3_core = x3[idx3]
y3_core = y3[idx3]

# Outlet section length (Core)
L_outlet_core = y3_core - y2_core

# Average Y at outlet (Core)
y_avg_core = (y2_core + y3_core) / 2

#Velocity at outlet (Scavenger and Core, we assumed that V_Scavenger = V_Core)
V_outlet = (V_inlet * L_inlet) / (L_outlet_scavenger + L_outlet_core)


# Mass flow rate at Scavenger outlet
m_dot_scavenger = rho * V_outlet * L_outlet_scavenger

# Mass flow rate at Core outlet
m_dot_core = rho * V_outlet * L_outlet_core




print(f"\nVelocity at Outlet (Scavenger) At x = {10} m, and y = {y_avg_scavenger:.3f} m: Vx = {V_outlet:.3f} m/s")

print(f"Mass flow rate at Scavenger outlet = {m_dot_scavenger:.4f} kg/s")

print(f"\nVelocity at Outlet (Core) At x = {10} m, and y = {y_avg_core:.3f} m: Vx = {V_outlet:.3f} m/s")

print(f"Mass flow rate at Core outlet = {m_dot_core:.4f} kg/s")
############################################################################################################################################
def generate_particles(num_particles=100):
    # Sizes and their percentages
    sizes = [1, 5, 10, 50, 100]  # micrometers
    percentages = [10, 30, 25, 20, 15]  # %

    particles = []

    # Generate based on percentages
    for size, perc in zip(sizes, percentages):
        count = int(num_particles * perc / 100)
        particles.extend([size] * count)

    # Shuffle so they're not grouped
    np.random.shuffle(particles)

    return np.array(particles)


distBetween = L_inlet/5
startPos = [x2_near,y2_near+distBetween]
max_offset = 0.6495 * 0.05
releasedParticles = generate_particles()

#for all sand grains released at a point run sim
for i in range(3):
    print(i)
    for j, pd in enumerate(releasedParticles):
        print(pd)
        print(startPos)
        # Random offset in range [-max_offset, +max_offset]
        rand_offset = np.random.uniform(-max_offset, max_offset)
        # Particle velocity = airflow + random offset in x only
        startVelocity = np.array([1.000,0.000]) + np.array([rand_offset, 0.0])
        #run sim
        trajectoryCalculations.run_particle_simulation(startPos,startVelocity,pd)
    #go up to next start position in y direction
    startPos = startPos[1] + distBetween

scavenge,core,lost = trajectoryCalculations.get_outlet_lists()
print(f"\nMasses in Scavenge", scavenge)
print(f"\nMasses in Core", core)
print(f"\nMasses Lost", lost)




