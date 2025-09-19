
import numpy as np
import MotionParticles

#time step for when we calculate each particles position
t_step = 0.01

etan = 0.70
enorm = 0.95



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


    return (abs(Vx), abs(Vy))


def update_boundary_matrix(matrix):
    boundary_matrix = matrix
    return

def run_particle_simulation(startPos,startVelocity,pd):
    #not sure what the return should be here yet

    particle_position_matrix =([startPos])
    pos = np.array(startPos)
    delta_time = 0
    collision = False
    Up = startVelocity
    counter = 0
    while (pos[0] <= (-0.3)): #run through sim until the particle crosses x boundary
        if (collision):
            Up = wall_relfection(Up)
            #print("collision vector", Up)
            #print("afvector", airflowVector)
            #print("dragVector", dragVector)
            collision = False
        airflowVector=np.array(MotionParticles.air_velocity(pos[0],pos[1]))

        dragVector = schiller_naumann_drag(Up,airflowVector,pd)

        Up = new_velocity(Up,dragVector)
        #print("sum of V",Up)
        pos = np.array(update_particle_position(pos, Up, t_step))
        #print("new pos", pos)
        particle_position_matrix = np.vstack([particle_position_matrix, pos])
        collision = collision_check(pos,pd)
        #print("collision", collision)
        delta_time+=t_step
    #if last pos is in scavange add to counter
    if (pos[1]<0):
        return True
    else:
        return False

def update_particle_position(startPos,Up,delta_time):
    #using new velocity calculate new position
    #Up is total particle velocity at this time
    #returns a 2-D array with the particles new [x,y]
    new_pos = startPos + (Up * delta_time)
    return np.array(new_pos)

def new_velocity(particleVector,dragVector):
    #calulate force balance
    Up = particleVector + dragVector * t_step
    return Up


import numpy as np


def schiller_naumann_drag(v, u, dp):
    #Compute drag acceleration on a 2D particle using Schiller-Naumann correlation.
    #Parameters:
    #v : particle velocity in microm/s
    #u  : local fluid (air) velocity in m/s
    #dp : particle diameter (m)
    #Returns: particle acceleration in m/s^2 : np.array([ax, ay])

    #particle density (kg/m^3)
    rho_p = 2650
    #fluid density (kg/m^3)
    rho_f = 1.184
    #dynamic viscosity of fluid (Pa·s)
    mu = 0.00001849

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


def has_collided(particlePos, pd):
    #checks if current particle position is in the boundary matrix
    #adds tolerance based on particle diameter (pd)
    #returns true if it has collided false if not
    #matches = np.all(boundary_matrix == particlePos, axis=1)
    matches = np.all(np.abs(boundary_matrix - particlePos) <= (pd*0.5), axis=1)
    return matches

def collision_check(particlePos, pd):
    dists = np.linalg.norm(boundary_matrix - particlePos, axis=1)
    nearest_idx = np.argmin(dists)
    min_dist = dists[nearest_idx]

    if min_dist <= pd / 2:
        return True
    else:
        return False

def wall_relfection(particle_vector):
    #switches the y when colliding with wall and takes into account given CR
    #recturns a vector [x,y]
    vxa = particle_vector[0] * enorm
    vya = particle_vector[1] * -etan
    particle_vector_after = np.array(vxa,vya)
    return particle_vector_after
