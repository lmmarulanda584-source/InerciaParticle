
import numpy as np
import math

def update_particle_position(startPos,Up):
    #using new velocity caluclate new positon
    new_pos = [0,0] #<-- this needs to be calculated
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
