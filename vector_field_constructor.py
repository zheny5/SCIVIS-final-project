import numpy as np
import sys

class vector_field():
    def __init__(self,dimensions=None, bounds=None):
        if dimensions is None:
            self._dimensions = [100, 100, 100]
        else:
            self._dimensions = dimensions
        if bounds is None:
            self._bounds = [-2.5, 2.5, -2.5, 2.5, 0., 8.]
        self.coords = list(np.linspace(self._bounds[2 * i], self._bounds[2 * i + 1], self._dimensions[i]) 
                      for i in range(3))
        self.grid = np.meshgrid(*self.coords, indexing='ij')

        self.u = self.skewing_oscillating_gyre_saddle(*self.grid)
    def skewing_oscillating_gyre_saddle(self, x, y, t):
        def gyre_saddle(a, x, y):
            def clamped_cos(x):
                return np.cos(np.clip(x, -np.pi / 2.0, np.pi / 2.0))

            def inner_saddle(x, y):
                dx = -np.sin(np.pi * x) * np.cos(np.pi * y) + a * np.sin(np.pi * y) * np.cos(np.pi * x)
                dy = np.sin(np.pi * y) * np.cos(np.pi * x) - a * np.sin(np.pi * x) * np.cos(np.pi * y)
                return dx, dy

            def outer_saddle(x, y):
                k = (y >= np.abs(x)).astype(float) - (y <= -np.abs(x)).astype(float)
                l = (x > np.abs(y)).astype(float) - (x < -np.abs(y)).astype(float)
                dx = k * a * clamped_cos(np.pi * x - a * np.pi * (y - k / 2)) \
                        - l * clamped_cos(np.pi * y - a * np.pi * (x - l / 2))
                dy = k * clamped_cos(np.pi * x - a * np.pi * (y - k / 2)) \
                        - l * a * clamped_cos(np.pi * y - a * np.pi * (x - l / 2))
                return dx, dy

            inside = (np.abs(x) <= 0.5) & (np.abs(y) <= 0.5)
            dx_inner, dy_inner = inner_saddle(x, y)
            dx_outer, dy_outer = outer_saddle(x, y)
            dx = np.where(inside, dx_inner, dx_outer)
            dy = np.where(inside, dy_inner, dy_outer)
            return dx, dy

        a = .5 * (np.sin(t * (2. * np.pi) / 1.0))

        x -= 0.2 * np.sin(t * (2.0 * np.pi) / 4.0)
        y -= 0.8 * np.sin(t * (2.0 * np.pi) / 4.0)

        dx, dy = gyre_saddle(a, x, y)
        dx *= 1.5
        dy *= 1.5
        dt = np.ones_like(t)

        return np.stack([dx, dy, dt], axis=-1)

def main():
    # construct analytical vector field
    vf = vector_field()
    print(vf.u.shape)

    

if __name__ == '__main__':
    main()
