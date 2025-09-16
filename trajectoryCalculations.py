
import numpy as np
import math

#time step for when we calculate each particles position
t_step = 0.01

boundary_matrix = np.array([[0, 0]])
def update_boundary_matrix(matrix):
    boundary_matrix = matrix
    print(boundary_matrix) #this print out is just for debugging

def run_particle_simulation(startPos,startVelocity, startAirFlow):
    #not sure what the return should be here yet
    particle_position_matrix =([startPos])
    pos = startPos
    delta_time = 0
    collision = False
    while (pos[0] <= 10): #run through sim until the particle crosses x boundary
        if (collision):
            wall_relfection() #not sure what's needed here Carolina to fill in this function
        Up = new_velocity(startVelocity,airflowVector=[0,3]) #<-- still need to properly calculate airflow vector
        pos = update_particle_position(startPos,Up,delta_time)
        particle_position_matrix = np.vstack([particle_position_matrix, pos])
        collision = has_collided(pos)
        delta_time+=t_step
    return particle_position_matrix

def update_particle_position(startPos,Up,delta_time):
    #using new velocity calculate new position
    #Up is total partical velocity at this time
    #returns a 2-D array with the particles new [x,y]
    #new_pos = startPos + (Up * delta_time) <-- todo: fix and seperate to x,y
    new_pos = [12,0]
    return new_pos

def new_velocity(startVector,airflowVector):
    #calulate force balance
    #still need to find the equation for this
    Up = 0  #<-- this needs to be calculated
    return Up

def air_velocity(theta,r):
    """
    Calculate air Velocity at any given point
    Parameter: theta,r
    Returns: Vr,Vtheta

    r: radial distance
    theta: angle in radians
    returns: (Vr, Vtheta)
    """
    Vr = -math.sin(theta)
    Vtheta = 1 / (2 * math.pi * r) + math.cos(theta)
    return Vr, Vtheta


def schiller_naumann_drag(u, up, rho_f, mu_f, dp):
    """
    Calculate drag force on a spherical particle using Schiller-Naumann model.

    Parameters
    ----------
    u : tuple (ux, uy)
        Fluid velocity vector
    up : tuple (upx, upy)
        Particle velocity vector
    rho_f : float
        Fluid density [kg/m^3]
    mu_f : float
        Fluid dynamic viscosity [Pa·s]
    dp : float
        Particle diameter [m]

    Returns
    -------
    Fx, Fy : float
        Drag force components [N]
    """
    # Relative velocity
    ur = (u[0] - up[0], u[1] - up[1])
    ur_mag = math.sqrt(ur[0] ** 2 + ur[1] ** 2)

    # Particle Reynolds number
    Re_p = rho_f * ur_mag * dp / mu_f

    # Schiller-Naumann drag coefficient
    if Re_p < 1000:
        Cd = 24 / Re_p * (1 + 0.15 * Re_p ** 0.687) if Re_p > 1e-8 else 0.0
    else:
        Cd = 0.44

    # Particle projected area (circle)
    Ap = math.pi * (dp ** 2) / 4

    # Drag force magnitude
    Fd_mag = 0.5 * Cd * rho_f * Ap * ur_mag ** 2

    # Force vector (same direction as relative velocity)
    if ur_mag > 1e-12:
        Fx = Fd_mag * ur[0] / ur_mag
        Fy = Fd_mag * ur[1] / ur_mag
    else:
        Fx, Fy = 0.0, 0.0

    return Fx, Fy

def has_collided(particlePos):
    #checks if currect particle positon is in the boundry matrix
    #returns true if it has collided false if not
    matches = np.all(boundary_matrix == particlePos, axis=1)
    return matches

def wall_relfection():
    #assumes it bounced to a point
    return
