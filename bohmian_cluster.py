import numpy as np
from particle import Particle
from functions import uniform_radial_sample


class BohmianCluster:
    # Sample particles in the initial probability distribution from a maximum radius
    def __init__(self, wavefunction, particle_count, sample_radius=10):
        self.particles = []
        self.wavefunction = wavefunction
        self.sample_radius = sample_radius

        while len(self.particles) < particle_count:
            batch = np.array([uniform_radial_sample(self.sample_radius) for _ in range(particle_count)]).T
            wavefunction.compute(batch[0], batch[1], batch[2], np.zeros(particle_count))
            prob_weight = 1 / max(self.wavefunction.density_value)
            self.particles += [Particle(batch[0][i], batch[1][i], batch[2][i]) for i in range(particle_count) if np.random.rand() / prob_weight < self.wavefunction.density_value[i]]
            print("samples so far: ", len(self.particles))

        if len(self.particles) > particle_count:
            self.particles = self.particles[:particle_count]

        print("sampling done!")

        self.xa = np.array([self.particles[i].x for i in range(len(self.particles))])
        self.ya = np.array([self.particles[i].y for i in range(len(self.particles))])
        self.za = np.array([self.particles[i].z for i in range(len(self.particles))])
        self.t = np.zeros(len(self.particles))


    # Slow but reliable update step for hydrogen atoms
    def rungekutta_step(self, timestep):
        self.wavefunction.compute(self.xa, self.ya, self.za, self.t)

        k1 = self.wavefunction.bohmian_trajectory * 1.0

        self.wavefunction.compute(self.xa + timestep * 0.5 * k1[0], self.ya + timestep * 0.5 * k1[1], self.za + timestep * 0.5 * k1[2], self.t + timestep * 0.5)
        k2 = self.wavefunction.bohmian_trajectory * 1.0

        self.wavefunction.compute(self.xa + timestep * 0.5 * k2[0], self.ya + timestep * 0.5 * k2[1], self.za + timestep * 0.5 * k2[2], self.t + timestep * 0.5)
        k3 = self.wavefunction.bohmian_trajectory * 1.0

        self.wavefunction.compute(self.xa + timestep * k3[0], self.ya + timestep * k3[1], self.za + timestep * k3[2], self.t + timestep)
        k4 = self.wavefunction.bohmian_trajectory * 1.0

        directions = timestep * (k1 + 2 * k2 + 3 * k3 + k4) / 6
        self.t += timestep
        return directions[0], directions[1], directions[2]
    

    # Fast but really bad update step for hydrogen atoms; orbitals of fixed/periodic rho are unstable with Euler updates
    def euler_step(self, timestep):
        self.wavefunction.compute(self.xa, self.ya, self.za, self.t)
        directions = timestep * self.wavefunction.bohmian_trajector

        self.t += timestep
        return directions[0], directions[1], directions[2]


    # Get the next timesteps and perform updates on the particles themselves.
    def update(self, timestep):
        dx, dy, dz = self.rungekutta_step(timestep)
        self.xa += np.real(dx)
        self.ya += np.real(dy)
        self.za += np.real(dz)
        
        for i in range(len(self.particles)):
            self.particles[i].update(self.xa[i], self.ya[i], self.za[i], True)
