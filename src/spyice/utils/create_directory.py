import os


def create_output_directory(
    hyd_dir: str,
    initial_salinity,
    boundary_condition_type,
    grid_resolution_dz,
    grid_timestep_dt,
    max_iterations,
    output_suffix,
) -> str:
    """Creates an output directory for storing temperature data.

    Args:
        hyd_dir (str): The directory path where the output directory will be created.
    Returns:
        str: The path of the created output directory.
    Raises:
        None
    """

    file_name = f"Temperature_{initial_salinity}_{boundary_condition_type}_{grid_resolution_dz}_{grid_timestep_dt}_{max_iterations}_{output_suffix}"
    output_dir = os.path.join(hyd_dir, file_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir
