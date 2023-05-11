import py_fast_marching_method as fmm
import numpy as np
import matplotlib.pyplot as plt


def signed_arrival_time_example():
    grid_size = np.array([50, 100])
    # grid_spacing = 1.0 / grid_size
    grid_spacing = np.ones(grid_size.shape)
    boundary_indices = np.array([[31, 75]])
    boundary_times = np.array([0.0])
    uniform_speed = 1.0

    arrival_times = fmm.uniform_speed_signed_arrival_time(
        grid_size, boundary_indices, boundary_times, grid_spacing, uniform_speed
    )

    print("Max dist:", np.max(arrival_times[:]))

    plt.imshow(arrival_times)
    plt.plot(boundary_indices[0,1], boundary_indices[0,0], "mo")
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    signed_arrival_time_example()
