import numpy as np

class RNGHandler:
    def __init__(self, seed=None):
        """
        Initialize the RNGHandler with an optional seed.
        """
        self.rng = np.random.default_rng(seed)

    def random(self):
        """Generate a random float in [0, 1)."""
        return self.rng.random()

    def uniform(self, low, high):
        """Generate a random float in the range [low, high)."""
        return self.rng.uniform(low, high)

    def log_uniform(self, scale):
        """Generate a random float for exponential distribution."""
        return -np.log(1 - self.rng.random()) / scale

    def choice(self, a, size=None, replace=True, p=None):
        """Generate random choices from a sequence."""
        return self.rng.choice(a, size=size, replace=replace, p=p)
