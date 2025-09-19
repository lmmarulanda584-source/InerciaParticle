import numpy as np
import matplotlib.pyplot as plt

def cot(x):
    return 1 / np.tan(x)

def air_velocity(x, y):
    # inlet at x = x_inlet
    if x <= -10:  # or within a small distance from inlet
        return 1.0, 0.0  # 1 m/s horizontal, no vertical component

    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    Vr = -np.sin(theta)
    Vt = (1 / (2 * np.pi * r)) + np.cos(theta)

    Vx = Vr * np.cos(theta) - Vt * np.sin(theta)
    Vy = Vr * np.sin(theta) + Vt * np.cos(theta)
    return Vx, Vy

def generate_curve(C, sign=1, n_points=5000):
    t = np.linspace(np.deg2rad(15.5568), np.deg2rad(176.891), n_points)
    if sign == 1:
        x = (5 * (C - t) * cot(t)) / (2*np.pi)
        y = 5 * (C - t) / (2*np.pi)
    else:
        x = (5 * (C - t) * cot(t)) / (2*np.pi)
        y = -5 * (C - t) / (2*np.pi)
    return np.column_stack((x, y))

# Curves
C1 = 3.7699111843077517
C2 = 4.39822971502571
curve1 = generate_curve(C1, sign=1)
curve2 = generate_curve(C2, sign=-1)
curve3 = generate_curve(np.pi, sign=-1)
curve4 = generate_curve(np.pi, sign=1)
walls = [curve1, curve2, curve3, curve4]

# Simulation parameters
dt = 0.1
steps = 30000
enormal = 0.95
etangential = 0.70

# Initial position
x, y = -8.0, 0.0
vx, vy = 1.0, 0.0  # initial velocity

trajectory_x = [x]
trajectory_y = [y]

for _ in range(steps):
    # Air velocity at the current position
    Vx_air, Vy_air = air_velocity(x, y)

    alpha = 0.1  # how fast the particle adjusts to the air
    vx += alpha * (Vx_air - vx)
    vy += alpha * (Vy_air - vy)

    # New position
    x_new = x + vx * dt
    y_new = y + vy * dt

    # Check collisions with walls
    for curve in walls:
        distances = np.sqrt((curve[:, 0] - x_new) ** 2 + (curve[:, 1] - y_new) ** 2)
        idx_min = np.argmin(distances)
        if distances[idx_min] < 0.1:  # collision threshold
            # Closest point
            px, py = curve[idx_min]

            # Tangent vector approximation (local derivative)
            if 1 < idx_min < len(curve) - 2:
                tx = curve[idx_min + 1, 0] - curve[idx_min - 1, 0]
                ty = curve[idx_min + 1, 1] - curve[idx_min - 1, 1]
            else:
                tx, ty = 1, 0
            tangent = np.array([tx, ty])
            tangent = tangent / np.linalg.norm(tangent)

            # Normal = perpendicular to the tangent
            normal = np.array([-tangent[1], tangent[0]])

            # Decompose velocity
            v = np.array([vx, vy])
            v_t = np.dot(v, tangent) * tangent
            v_n = np.dot(v, normal) * normal

            # Apply coefficients
            v_rebound = etangential * v_t - enormal * v_n
            vx, vy = v_rebound

            # Push slightly outside the wall
            x_new = x + vx * dt + normal[0] * 1e-4
            y_new = y + vy * dt + normal[1] * 1e-4

            break

    # Save trajectory
    trajectory_x.append(x_new)
    trajectory_y.append(y_new)

    x, y = x_new, y_new

# ---------------------------
# Plot
# ---------------------------
plt.figure(figsize=(8,6))
plt.plot(curve1[:,0], curve1[:,1], 'k')
plt.plot(curve2[:,0], curve2[:,1], 'k')
plt.plot(curve3[:,0], curve3[:,1], 'k')
plt.plot(curve4[:,0], curve4[:,1], 'k')
plt.scatter(trajectory_x, trajectory_y, c='r', s=10, label="Trajectory")

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Trajectory of a particle with rebounds")
plt.legend()
plt.axis("scaled")
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.show()





