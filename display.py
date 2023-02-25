import numpy as np
from tkinter import *


# Take a point from the system's 3d space and project it onto the 2d space to render
def project(x, y, z, dim, t=0, period=1):
    x = x * np.cos(t / period) + z * np.sin(t / period)
    z = x * -np.sin(t / period) + z * np.cos(t / period)
    scale = 0.15 / max(0.2, z * 0.2 + 5)
    newX = x * scale * dim + dim / 2
    newY = y * scale * dim + dim / 2
    newZ = z
    return newX, newY, newZ


# Add a circle function to the canvas
def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x-r, y-r, x+r, y+r, **kwargs)
Canvas.create_circle = _create_circle


# A container class to handle a 3d cluster and a 2d TKinter display
class Display:
    def __init__(self, cluster, root, canvas, height, width, fps=60, lifetime=1000):
        self.root = root
        self.canvas = canvas
        self.height = height
        self.width = width
        self.realtime_tick = int(1000 / fps)
        self.cluster = cluster
        self.simulation_tick = 0.2

        self.lifetime = lifetime
        self.timer = 25
        self.rotation_period = 1e9
        
    # Single frame of display render
    def render(self):
        if self.timer < self.lifetime:
            c = self.canvas
            c.delete("all")

            # Render reference points (center) and information
            c.create_text(100,10,fill="white",font="Arial 12",
                        text="simulation progress: " + str(self.timer) + "/" + str(self.lifetime))
            cx, cy, cz = project(0, 0, 0, self.height, self.timer, self.rotation_period)
            c.create_circle(cx, cy, 2, fill="red", outline="black", width=0)

            # Each particle has a buffer of points that we'll only render the first of for now.
            for p in self.cluster.particles:
                x0, y0, z0 = project(p.buffer[0][0], p.buffer[0][1], p.buffer[0][2], self.height, self.timer, self.rotation_period)
                c.create_circle(x0, y0, 1, fill="green", outline="black", width=0)

            # Update the state of the atom
            self.cluster.update(self.simulation_tick)
            self.timer += 1
                
            # Callback loop
            self.canvas.after(self.realtime_tick, self.render)


