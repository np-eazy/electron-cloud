from wavefunction import *
from bohmian_cluster import *
from tkinter import *
from display import Display

# The eigenstates governing this instance of a hydrogen atom
wavefunction = Wavefunction([Eigenstate(1, 0, 0, 0.707), Eigenstate(2, 1, 1, 0.707j)])

# The number of particles to simulate. 
cluster = BohmianCluster(wavefunction, 100)

# Display dimensions
DIM = 500
MAX_FRAMERATE = 60
FRAMES = 2000



# Execution starts
xDim, yDim = DIM, DIM
root = Tk()
canvas = Canvas(root, width=xDim, height=yDim, borderwidth=0, highlightthickness=0,
                   bg="black")
canvas.grid()
''
display = Display(cluster, root, canvas, xDim, yDim, MAX_FRAMERATE, FRAMES)
display.render()
root.title("Hydrogen Atom")
root.mainloop()