from dataclasses import dataclass
import numpy as np


def set_dataclass(data_to_be_converted: dict, dataclass: dataclass):
    """Sets the values of a dataclass object using a dictionary.

    Args:
        data_to_be_converted (dict): The dictionary containing the data to be converted.
        dataclass (dataclass): The dataclass object to be updated.
    Returns:
        dataclass: The updated dataclass object.
    """

    for key, value in data_to_be_converted.items():
        setattr(dataclass, key, value)

    return dataclass


### Computes total time passed
def t_total(t_passed, dt):
    """Computes total time passed based on current time step dt and, total time of previous time step.

    Args:
        t_passed (float): The time that has already passed.
        dt (float): The time increment.
    Returns:
        float: The total time.
    """

    return t_passed + dt


def export_residuals(
    self, residuals: list, temperature_mushy, phi_mushy, salinity_mushy, output_dir
):
    """Exports the residuals to a file.

    Args:
        residuals (np.array): The residuals to export.
        output_dir (str): The output directory to save the residuals to.
    Returns:
        None
    """

    residuals = np.array(residuals, dtype=object)
    temperature_mushy = np.array(temperature_mushy, dtype=object)
    phi_mushy = np.array(phi_mushy, dtype=object)
    salinity_mushy = np.array(salinity_mushy, dtype=object)
    np.save(output_dir + "/residuals.npy", residuals)
    np.save(output_dir + "/temperature_mushy.npy", temperature_mushy)
    np.save(output_dir + "/phi_mushy.npy", phi_mushy)
    np.save(output_dir + "/salinity_mushy.npy", salinity_mushy)
    print("Residuals exported successfully.")


### Initialize time step, starting time, end time


# def CFL_time_step(w,dt,dz,nz):
#     '''
#     time step advection
#     '''
#     dt = np.zeros(nz)
#     for i in range(nz):
#         dt[i] = abs(dz/w[i])
#     dt_non_zero = dt[dt !=0]
#     dt_CFL = np.amin(dt_non_zero)
#     return dt_CFL
