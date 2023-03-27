import py_fast_marching_method as fmm
import numpy as np
import matplotlib.pyplot as plt


def signed_arrival_time_example():
    grid_size = np.array([5, 5])
    grid_spacing = 1.0 / grid_size
    boundary_indices = np.array([[2, 2]])
    boundary_times = np.array([0.0])
    uniform_speed = 1.0

    arrival_times = fmm.uniform_speed_eikonal_signed_arrival_time(
        grid_size, boundary_indices, boundary_times, grid_spacing, uniform_speed
    )
    plt.imshow(arrival_times)
    plt.show()


if __name__ == "__main__":
    signed_arrival_time_example()
