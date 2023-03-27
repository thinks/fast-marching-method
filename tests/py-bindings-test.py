import py_fast_marching_method as fmm
import numpy as np
from scipy import ndimage as ndi


def test_main():
    grid_size = np.array([5, 5])
    # grid_spacing = 1.0/grid_size
    grid_spacing = np.ones(grid_size.shape)
    boundary_indices = np.array([[2, 2]])
    boundary_times = np.array([0.0])
    uniform_speed = 1.0

    arrival_times = fmm.uniform_speed_eikonal_signed_arrival_time(
        grid_size, boundary_indices, boundary_times, grid_spacing, uniform_speed
    )
    print("arrival_times\n", arrival_times)

    binary_input_grid = np.ones(grid_size)
    binary_input_grid[boundary_indices[:, 0], boundary_indices[:, 1]] = 0
    print("binary_input_grid\n", binary_input_grid)

    edt = ndi.distance_transform_edt(binary_input_grid)
    print("edt\n", edt)

    print("All clear")


if __name__ == "__main__":
    test_main()
