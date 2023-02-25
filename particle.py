# A basic particle container class containing a buffer in case we want to render trails
class Particle:
    def __init__(self, x, y, z, maxBufLength=6):
        self.x = x
        self.y = y
        self.z = z
        self.buffer = [[x, y, z] for _ in range(maxBufLength)]
        self.maxBufLength = maxBufLength

    # If set_values is False, we treat x y z as increments to change the currenct values.
    def update(self, x, y, z, set_values=False):
        if set_values:
            self.x = x
            self.y = y
            self.z = z
        else:
            self.x += x
            self.y += y
            self.z += z
        self.buffer = self.buffer[1:] + [[self.x, self.y, self.z]]