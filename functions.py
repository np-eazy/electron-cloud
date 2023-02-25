import numpy as np
import scipy

# Special functions
def associated_laguerre(n, l, rad):
    norm = np.sqrt(np.power(2 / (n * rad), 3) * scipy.special.factorial(n - l - 1) / (2 * n * scipy.special.factorial(n + l)))
    def func(rho):
        degree = n - l - 1
        alpha = 2 * l + 1
        dropoff = 2 * rho / (n * rad)
        return scipy.special.assoc_laguerre(rho, degree, alpha) * norm * np.exp(dropoff / -2) * np.power(rho, l)
    return func

def associated_legendre(l, m):
    def func(theta):
        if theta.shape == ():
            return scipy.special.lpmn(m, l, np.cos(theta))[0][-1, -1]
        else:
            return np.array([scipy.special.lpmn(m, l, np.cos(t))[0][-1, -1] for t in theta])
    return func

def polar_harmonic(m):
    def func(phi):
        return np.exp(1j * m * phi)
    return func
    
# Spherical utils
def secant_line(x, func, h=1e-8):
    return (func(x + h/2) - func(x - h/2)) / h

def spherical_map(x, y, z):
    rho = np.sqrt(x * x + y * y + z * z)
    plane = np.sqrt(x * x + y * y)
    sign = y / np.abs(y)
    theta = np.arctan(plane / z)
    phi = np.arccos(x / plane) * sign
    return np.array([rho, theta, phi])

def spherical_basis(rho, theta, phi):
    rho_basis = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    theta_basis = np.array([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)])
    phi_basis = np.array([-np.sin(phi), np.cos(phi), 0 * phi])
    
    return np.array([rho_basis, theta_basis, phi_basis])


def spherical_gradient(rho, theta, phi, rho_func, theta_func, phi_func, rho_eval, theta_eval, phi_eval, h=1e-9):

    rho_partial = secant_line(rho, rho_func, h) * theta_eval * phi_eval
    theta_partial = secant_line(theta, theta_func, h) * phi_eval * rho_eval
    phi_partial = secant_line(phi, phi_func, h) * rho_eval * theta_eval

    coefs = np.array([rho_partial, theta_partial / rho, phi_partial / (rho * np.sin(theta))])
    basis = spherical_basis(rho, theta, phi)
    rho_unit, theta_unit, phi_unit = basis[0], basis[1], basis[2]
    
    return rho_unit * coefs[0] + theta_unit * coefs[1] + phi_unit * coefs[2]


def uniform_radial_sample(max_rad=5):
    rho, phi = np.sqrt(np.random.rand()) * max_rad, np.random.rand() * 2 * np.pi
    theta = np.arccos(1 - 2 * np.random.rand())
    return rho * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])




    









        


