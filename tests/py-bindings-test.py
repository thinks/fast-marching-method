import py_fast_marching_method as fmm
import numpy as np
from scipy import ndimage as ndi


def test_uniform_2D():
    grid_size = np.array([101, 101])
    # grid_spacing = 1.0/grid_size
    grid_spacing = np.ones(grid_size.shape)
    boundary_indices = np.array([[50, 50]])
    boundary_times = np.array([0.0])
    uniform_speed = 1.0

    arrival_times = fmm.uniform_speed_signed_arrival_time(
        grid_size, boundary_indices, boundary_times, grid_spacing, uniform_speed
    )
    print("arrival_times\n", arrival_times)

    binary_input_grid = np.ones(grid_size)
    binary_input_grid[boundary_indices[:, 0], boundary_indices[:, 1]] = 0
    print("binary_input_grid\n", binary_input_grid)

    edt = ndi.distance_transform_edt(binary_input_grid)
    print("edt\n", edt)

    # FMM is not super accurate -> large tolerance
    np.testing.assert_allclose(arrival_times, edt, rtol=0, atol=0.5)

    print("Done 2D")


def test_uniform_3D():
    grid_size = np.array([101, 101, 101])
    # grid_spacing = 1.0/grid_size
    grid_spacing = np.ones(grid_size.shape)
    boundary_indices = np.array([[50, 50, 50]])
    boundary_times = np.array([0.0])
    uniform_speed = 1.0

    arrival_times = fmm.uniform_speed_signed_arrival_time(
        grid_size, boundary_indices, boundary_times, grid_spacing, uniform_speed
    )
    print("arrival_times\n", arrival_times)

    binary_input_grid = np.ones(grid_size)
    binary_input_grid[
        boundary_indices[:, 0], boundary_indices[:, 1], boundary_indices[:, 2]
    ] = 0
    print("binary_input_grid\n", binary_input_grid)

    edt = ndi.distance_transform_edt(binary_input_grid)
    print("edt\n", edt)

    # FMM is not super accurate -> large tolerance
    np.testing.assert_allclose(arrival_times, edt, rtol=0, atol=0.7)

    print("Done 3D")


def test_varying_2D():
    grid_size = np.array([5, 5])
    # grid_spacing = 1.0/grid_size
    grid_spacing = np.ones(grid_size.shape)
    boundary_indices = np.array([[2, 2]])
    boundary_times = np.array([0.0])
    varying_speed = np.ones(grid_size)

    arrival_times = fmm.varying_speed_signed_arrival_time(
        grid_size, boundary_indices, boundary_times, grid_spacing, varying_speed
    )
    print("arrival_times\n", arrival_times)

    binary_input_grid = np.ones(grid_size)
    binary_input_grid[boundary_indices[:, 0], boundary_indices[:, 1]] = 0
    print("binary_input_grid\n", binary_input_grid)

    edt = ndi.distance_transform_edt(binary_input_grid)
    print("edt\n", edt)

    # FMM is not super accurate -> large tolerance
    np.testing.assert_allclose(arrival_times, edt, rtol=0, atol=0.5)

    print("Done 2D")


if __name__ == "__main__":
    test_uniform_2D()
    test_uniform_3D()
    test_varying_2D()
