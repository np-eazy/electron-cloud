Atom Visualization using Bohmian Trajectories
by Joey Zhu (github username @np-eazy)

2/24/23: Right now there is a fully functioning Python implementation that shows the simulations for each eigenstate are accurate, but the rendering in tkinter has very rudimentary and slow functions, and the computation can be better organized. I have no intention of taking the Python implementation any further; in the future, I plan to write a memory-level implementation based on 'Parallelization Blueprint.xlsx' to better improve the performance. It would also be much better to export the 3d points to a 3rd-party render engine like Blender.

Run: in the electron-cloud folder, run '$ python3 main.py' to generate a TKinter window that will start a new simulation.

See atom-2.mp4 for an example of what a render might look like; I forgot the magnitudes/phases but that one is created with (1, 0, 0), (2, 0, 0), (3, 0, 0), and (5, 1, -1) if I remember correctly. I changed the initialization code in the BohmianCluster class to get that distribution to start from a small point rather than to sample everything.

The projection code also pre-rotates along the y axis but this can definitely be customized.
