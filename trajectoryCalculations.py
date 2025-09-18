
import numpy as np

#time step for when we calculate each particles position
t_step = 0.01

etan = 0.70
enorm = 0.95

#counters
scavenge =np.array([])
core = np.array([])
lost = np.array([])

boundary_matrix = np.array([[0, 0]])

def air_velocity(x,y):
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

def get_outlet_lists():
    return scavenge,core,lost

def update_boundary_matrix(matrix):
    boundary_matrix = matrix
    return

def run_particle_simulation(startPos,startVelocity,pd):
    #not sure what the return should be here yet

    particle_position_matrix =([startPos])
    pos = startPos
    delta_time = 0
    collision = False
    Up = startVelocity
    while (pos[0] <= 10): #run through sim until the particle crosses x boundary
        if (collision):
            wall_relfection(Up)
            collision = False
        airflowVector=np.array(air_velocity(startPos[0],startPos[1]))
        dragVector = schiller_naumann_drag(Up,airflowVector,pd)
        Up = new_velocity(Up,dragVector)
        pos = update_particle_position(startPos,Up,delta_time)
        particle_position_matrix = np.vstack([particle_position_matrix, pos])
        collision = has_collided(pos)
        delta_time+=t_step
    #if last pos is in scavange add to counter
    if (pos[0]>2.3):
        np.append(scavenge,pd)
        #if core add to counter
    if (pos[0]<(-2.4)):
        #in neither add to counter
        np.append(core,pd)
    else:
        np.append(lost,pd)

    return

def update_particle_position(startPos,Up,delta_time):
    #using new velocity calculate new position
    #Up is total particle velocity at this time
    #returns a 2-D array with the particles new [x,y]
    new_pos = startPos + (Up * delta_time)
    return new_pos

def new_velocity(particleVector,dragVector):
    #calulate force balance
    Up = particleVector + dragVector * t_step
    return Up


import numpy as np


def schiller_naumann_drag(v, u, dp):
    """
    Compute drag acceleration on a 2D particle using Schiller-Naumann correlation.

    Parameters:
    v      : np.array([vx, vy]), particle velocity in microm/s
    u      : np.array([ux, uy]), local fluid (air) velocity in m/s
    dp     : float, particle diameter (m)

    Returns:
    a      : np.array([ax, ay]), particle acceleration in m/s^2
    """

    #particle density (kg/m^3)
    rho_p = 2650
    #fluid density (kg/m^3)
    rho_f = 1.225
    #dynamic viscosity of fluid (Pa·s)
    mu = .00001789

    #convert micro m to m
    dp = dp * 1e-6

    # relative velocity
    v_rel = u - v
    v_rel_mag = np.linalg.norm(v_rel)

    # avoid division by zero
    if v_rel_mag < 1e-12:
        return np.zeros_like(v)

    # Reynolds number
    Re_p = (rho_f * dp * v_rel_mag) / mu

    # Schiller–Naumann correction factor
    if Re_p < 1000:
        f = 1.0 + 0.15 * (Re_p ** 0.687)
    else:
        # fallback: constant drag coefficient ~0.44 for turbulent regime
        f = 0.44 * Re_p / 24.0

    # drag acceleration (per unit mass)
    a = (18.0 * mu / (rho_p * dp ** 2)) * f * v_rel
    return a


def has_collided(particlePos):
    #checks if current particle position is in the boundary matrix
    #returns true if it has collided false if not
    matches = np.all(boundary_matrix == particlePos, axis=1)
    return matches

def wall_relfection(particle_vector):
    vxa = particle_vector[0] * etan
    vya = particle_vector[1] * -enorm
    particle_vector_after = np.array(vxa,vya)
    return particle_vector_after
