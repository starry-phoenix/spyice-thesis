from dataclasses import dataclass


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
