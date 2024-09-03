import pytest

from spyice.parameters import UserInput


@pytest.mark.parametrize(
    "grid_time_step, initial_salinity, max_iterations",
    [(10, 3, 100), (30, 100, 5), (20, 3, 1000), (1, 100, 5)],
)
class TestUI:
    def test_userinput(self, grid_time_step, initial_salinity, max_iterations) -> None:
        ui = UserInput(
            grid_timestep_dt=grid_time_step,
            boundary_salinity=initial_salinity,
            max_iterations=max_iterations,
        )
        assert (
            ui.max_iterations > 1 & ui.max_iterations < 30000
        ), "iter_max must be greater than 0"
        assert ui.grid_timestep_dt > 0, "dt must be greater than 0"
        assert (
            ui.boundary_salinity > 1 & initial_salinity < 250
        ), "S_IC must be greater than 0 and less than 250"
