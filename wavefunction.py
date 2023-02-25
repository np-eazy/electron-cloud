from functions import *

# Quantum physical constants
HBAR = 2
ELECTRON_MASS = 1
ELECTRON_CHARGE = -1
PERMITTIVITY = 0.1
BOHR_RADIUS = 1
RYDBERG_CONST = ELECTRON_MASS * (ELECTRON_CHARGE ** 4) / (2 * ((4 * np.pi * PERMITTIVITY * HBAR) ** 2))

# Eigenstates and Wavefunctions can take in coordinates grouped by dimension instead by particle for vectorized computations.
class Eigenstate:
    def __init__(self, n, l, m, coef=1):
        self.n = n
        self.l = l
        self.m = m
        self.coef = coef

        # Functions associated with this eigenstate; associated_laguerre actually combines other radial functions
        self.rho_func = associated_laguerre(n, l, BOHR_RADIUS)
        self.theta_func = associated_legendre(l, m)
        self.phi_func = polar_harmonic(m)

        # The computational structure highly relies on recycling computations in each coordinate
        self.rho_value = None
        self.theta_value = None
        self.phi_value = None
        self.phase_value = 0
        self.wave_value = 0
        self.density_value = 0
        self.gradient_value = None
        self.antigradient = None

        self.flags = {
            "compute-wave": False,
            "compute-probability": False,
            "compute-gradient": False,
        }

    # Use the formula for a Hydrogen(ic) atom
    def get_energy(self, n):
        return -RYDBERG_CONST / (n * n)

    # Calculate the time-dependent phase of this eigenstate
    def get_phase(self, t):
        self.phase_value = self.coef * np.exp(1j * self.get_energy(self.n) * t / HBAR)
        return self.phase_value

    # Compute the wavefunction and set up to compute for probability densities and gradient
    def compute_wave(self, x, y, z, t):
        mapping = spherical_map(x, y, z)
        self.rho_value, self.theta_value, self.phi_value = mapping[0], mapping[1], mapping[2]
        self.rho_eval = self.rho_func(self.rho_value)
        self.theta_eval = self.theta_func(self.theta_value)
        self.phi_eval = self.phi_func(self.phi_value)

        self.wave_value = self.rho_eval * self.theta_eval * self.phi_eval * self.get_phase(t)
        self.flags["compute-wave"] = True
        self.flags["compute-density"] = False
        self.flags["compute-gradient"] = False
        return self.wave_value

    def compute_density(self):
        assert self.flags["compute-wave"]
        self.density_value = self.wave_value * np.conjugate(self.wave_value)
        self.flags["compute-density"] = True
        return self.density_value

    def compute_gradient(self):
        assert self.flags["compute-wave"]
        self.gradient_value = self.phase_value * spherical_gradient(self.rho_value, self.theta_value, self.phi_value, self.rho_func, self.theta_func, self.phi_func, self.rho_eval, self.theta_eval, self.phi_eval)
        self.flags["compute-gradient"] = True
        return self.gradient_value




class Wavefunction:
    def __init__(self, eigenstates):
        self.eigenstates = eigenstates
        self.wave_value = 0
        self.density_value = 0
        self.gradient_value = None

    def compute(self, x, y, z, t):
        # Aggregate superpositions of waves and gradients and use total wavefunction to compute probability density
        self.wave_value = np.sum(np.array([eigenstate.compute_wave(x, y, z, t) for eigenstate in self.eigenstates]), axis=0)
        self.gradient_value = np.sum(np.array([eigenstate.compute_gradient() for eigenstate in self.eigenstates]), axis=0)
        self.density_value = (self.wave_value * np.conjugate(self.wave_value)).real

        # Calculate the wavefunction
        self.bohmian_trajectory = ((np.conjugate(self.wave_value) * self.gradient_value) - (np.conjugate(self.gradient_value) * self.wave_value)) / self.density_value
        self.bohmian_trajectory = (self.bohmian_trajectory * (HBAR / (2j * ELECTRON_MASS))).real

