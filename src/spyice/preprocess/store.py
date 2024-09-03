from __future__ import annotations

import datetime
import os
import sys

import hickle as hkl
import numpy as np


def set_up_matrices(iter_max, nz):
    """Initialize matrices with one column for each time step

    Args:
        iter_max (int): The maximum number of iterations
        nz (int): The number of columns in the matrices
    Returns:
        tuple: A tuple containing the initialized matrices:
            - all_T (ndarray): A 2D array of shape (iter_max, nz) filled with zeros
            - all_S_sw (ndarray): A 2D array of shape (iter_max, nz) filled with zeros
            - all_phi (ndarray): A 2D array of shape (iter_max, nz) filled with zeros
            - all_H (ndarray): A 2D array of shape (iter_max, nz) filled with zeros
            - all_H_solid (ndarray): A 2D array of shape (iter_max, nz) filled with zeros
            - all_w (ndarray): A 2D array of shape (iter_max, nz) filled with zeros
            - all_thick (ndarray): A 1D array of length iter_max filled with zeros
            - all_t_passed (ndarray): A 1D array of length iter_max filled with zeros
    """
    all_T = np.zeros([iter_max, nz])
    all_S_sw = np.zeros([iter_max, nz])
    all_phi = np.zeros([iter_max, nz])
    all_H = np.zeros([iter_max, nz])
    all_H_solid = np.zeros([iter_max, nz])
    all_w = np.zeros([iter_max, nz])
    all_thick = np.zeros([iter_max])
    all_t_passed = np.zeros(iter_max)
    return all_T, all_S_sw, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed


def store_results(
    T,
    S_sw,
    phi,
    H,
    H_solid,
    w,
    thickness,
    t_passed,
    all_T,
    all_S_sw,
    all_phi,
    all_H,
    all_H_solid,
    all_w,
    all_thick,
    all_t_passed,
    t,
):
    """Stores the results of the simulation.
    Args:
        T (numpy.ndarray): Temperature values.
        S_sw (numpy.ndarray): S_sw values.
        phi (numpy.ndarray): phi values.
        H (numpy.ndarray): H values.
        H_solid (numpy.ndarray): H_solid values.
        w (numpy.ndarray): w values.
        thickness (float): Thickness value.
        t_passed (float): t_passed value.
        all_T (numpy.ndarray): Array to store all T values.
        all_S_sw (numpy.ndarray): Array to store all S_sw values.
        all_phi (numpy.ndarray): Array to store all phi values.
        all_H (numpy.ndarray): Array to store all H values.
        all_H_solid (numpy.ndarray): Array to store all H_solid values.
        all_w (numpy.ndarray): Array to store all w values.
        all_thick (numpy.ndarray): Array to store all thickness values.
        all_t_passed (numpy.ndarray): Array to store all t_passed values.
        t (int): Current time step.
    Returns:
        Tuple[numpy.ndarray]: A tuple containing all the updated arrays.
    """

    all_T[t, :] = T
    all_S_sw[t, :] = S_sw
    all_phi[t, :] = phi
    all_H[t, :] = H
    all_H_solid[t, :] = H_solid
    all_w[t, :] = w
    all_thick[t] = thickness
    all_t_passed[t] = t_passed
    return all_T, all_S_sw, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed


def create_directory():
    """Creates a new directory with the current timestamp as the directory name.

    Args:
        None
    Returns:
        None
    Raises:
        None
    """

    # Parent Directory path
    parent_dir = "/Users/asimson/work/GIT_projects/seaicemodel/data"  # "C:\\Users\\Anna Lara\\Documents\\05_GIT\\Seaicemodel\\Data\\"
    # Directory

    directory = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    # Path
    path = os.path.join(parent_dir, directory)
    os.mkdir(path)
    print("Directory '% s' created" % directory)
    sys.path.append(path)


def save_results(
    all_T, all_S, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed
):
    """Save results as hkl files.

    Args:
        all_T (list): List of temperatures.
        all_S (list): List of salinities.
        all_phi (list): List of liquid fractions.
        all_H (list): List of enthalpies.
        all_H_solid (list): List of solid state enthalpies.
        all_w (list): List of velocities.
        all_thick (list): List of thicknesses.
        all_t_passed (list): List of time passed.
    Returns:
        None
    """
    # go to path
    os.chdir(sys.path[-1])
    # save data
    hkl.dump(all_T, "temperature.hkl", compression="gzip")
    print("saved temperature")
    hkl.dump(all_phi, "liquid_fraction.hkl", compression="gzip")
    print("saved liquid fraction")
    hkl.dump(all_S, "salinity.hkl", compression="gzip")
    print("saved salinity")
    hkl.dump(all_H, "enthalpy.hkl", compression="gzip")
    print("saved enthalpy")
    hkl.dump(all_H_solid, "enthalpy_solid.hkl", compression="gzip")
    print("saved solid state enthalpy")
    hkl.dump(all_w, "velocity.hkl", compression="gzip")
    print("saved velocity")
    hkl.dump(all_thick, "thickness.hkl", compression="gzip")
    print("saved thickness")
    hkl.dump(all_t_passed, "time_passed.hkl", compression="gzip")
    print("saved time passed")
