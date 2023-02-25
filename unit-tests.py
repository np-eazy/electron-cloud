import functions
import numpy as np

def test_basis_functions():
    x = np.linspace(-2, 2, 500)
    for i in x:
        print(functions.associated_legendre(2, 0)(i))
    #plt.plot(x, functions.associated_legendre(2, 0))
    #plt.show()
    return None

def test_eigenstate():
    eigenstate = functions.Eigenstate(1, 0, 0, 1)
    print("eigenstate 1,0,0 at 1, 1, 1, 1: ", eigenstate.compute_wave(1, 1, 1, 1))
    print("probability density: ", eigenstate.compute_density())
    print("gradient: ", eigenstate.compute_gradient())
    print("Eigenstates passed")
    return None

def test_wavefunction():
    groundState = functions.Eigenstate(1, 0, 0, 0.707)
    excitedState = functions.Eigenstate(3, 2, 2, 0.707)
    wavefunction = functions.Wavefunction([groundState, excitedState])
    wavefunction.compute(1, 1, 1, 1)
    print("wavefunction at 1, 1, 1, 1: ", wavefunction.wave_value)
    print("probability density: ", wavefunction.density_value)
    print("gradient: ", wavefunction.gradient_value)
    print("Wavefunction passed")

def test_cluster():
    groundState = functions.Eigenstate(1, 0, 0, 0.707)
    excitedState = functions.Eigenstate(3, 2, 2, 0.707)
    wavefunction = functions.Wavefunction([groundState, excitedState])
    cluster = functions.BohmianCluster(wavefunction, 10)
    print("cluster particles: ", cluster.particles)
    for i in range(10):
        cluster.update(1)
        if i % 1 == 0:
            print(cluster.particles[0].x)
    print("Cluster passed")

test_basis_functions()
test_eigenstate()
test_wavefunction()
test_cluster()
