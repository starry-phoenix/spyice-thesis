from __future__ import annotations

from ..parameters import UserInput
from .initial_boundary_conditions import temperature_gradient

ui = UserInput()
(
    temperature_melt,
    boundary_top_temperature,
    boundary_salinity,
    boundary_top_temperature,
    dz,
    phi_c,
) = (
    ui.temperature_melt,
    ui.boundary_top_temperature,
    ui.boundary_salinity,
    ui.boundary_top_temperature,
    ui.grid_resolution_dz,
    ui.critical_liquid_fraction,
)


class ModifyInitialBoundary:
    """Defines functions for managing initial boundary conditions in a sea ice model."""

    def __init__(self) -> None:
        pass

    def bc_neumann(self, phi_k, nz, bc_condition=None):
        """Applies Neumann boundary condition to the given field.

        Args:
            phi_k (numpy.ndarray): The field to which the boundary condition is applied.
            nz (int): The number of grid points in the vertical direction.
            bc_condition (str, optional): The type of boundary condition. Defaults to None.
        Returns:
            None
        """

        if bc_condition == "Neumann":
            self.temp_grad = temperature_gradient(
                phi_k, nz
            )  # T_bottom - T_top / depth according to Buffo
        else:
            self.temp_grad = None

    def set_boundary_condition_type(self, critical_depth, bc_type="Neumann"):
        """Sets the boundary condition type for the given critical depth.

        Args:
            critical_depth (float): The critical depth value.
            bc_type (str, optional): The type of boundary condition. Defaults to "Neumann".
        Raises:
            None
        Returns:
            None
        """

        if bc_type == "Neumann":
            self.temp_grad = (
                dz * (temperature_melt - boundary_top_temperature) / critical_depth
            )
        else:
            self.temp_grad = None
