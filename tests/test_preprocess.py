import pytest

from spyice.parameters import Constants
from spyice.parameters.user_input import UserInput
from spyice.preprocess.pre_process import PreProcess


@pytest.mark.parametrize(
    "config_data, output_dir, expected_constants_type, expected_time_step, expected_initial_salinity, expected_max_iterations",
    [
        (
            {
                "constants_type": "real",
                "time_step": 1,
                "initial_salinity": 35,
                "max_iterations": 100,
            },
            "output/dir",
            Constants.REAL.value,
            1,
            35,
            100,
        ),
        (
            {"constants_type": "debug", "time_step": 2},
            "output/dir",
            Constants.DEBUG.value,
            2,
            None,
            None,
        ),
    ],
    ids=["real_constants", "debug_constants"],
)
def test_set_preprocess(
    config_data,
    output_dir,
    expected_constants_type,
    expected_time_step,
    expected_initial_salinity,
    expected_max_iterations,
):
    # Act
    preprocess = PreProcess()
    preprocess.set_preprocess(config_data, output_dir)

    # Assert
    assert preprocess.constants == expected_constants_type
    assert preprocess.dt == expected_time_step
    if expected_initial_salinity is not None:
        assert preprocess.S_IC == expected_initial_salinity
    if expected_max_iterations is not None:
        assert preprocess.iter_max == expected_max_iterations


@pytest.mark.parametrize(
    "config_data, output_dir",
    [
        (
            {
                "constants_type": "real",
                "time_step": 1,
                "initial_salinity": 35,
                "max_iterations": 100,
            },
            "output/dir",
        ),
        ({"constants_type": "debug", "time_step": 2}, "output/dir"),
    ],
    ids=["real_constants", "debug_constants"],
)
def test_preprocess(config_data, output_dir):
    # Arrange
    preprocess = PreProcess()
    preprocess.set_preprocess(config_data, output_dir)

    # Act
    preprocess.preprocess()

    # Assert
    assert hasattr(preprocess, "time_passed")
    assert hasattr(preprocess, "temperature")
    assert hasattr(preprocess, "salinity")
    assert hasattr(preprocess, "liquid_fraction")
    assert hasattr(preprocess, "omega")
    assert hasattr(preprocess, "solid_enthalpy")
    assert hasattr(preprocess, "enthalpy")
    assert hasattr(preprocess, "ice_thickness")


@pytest.mark.parametrize(
    "config, out_dir_final, expected_var1",
    [
        (
            {
                "constants_type": "real",
                "time_step": 1,
                "initial_salinity": 35,
                "max_iterations": 100,
            },
            "output/dir",
            "expected_value1",
        ),
        ({"constants_type": "debug", "time_step": 2}, "output/dir", "expected_value2"),
    ],
    ids=["real_constants", "debug_constants"],
)
def test_get_variables(config, out_dir_final, expected_var1):
    # Arrange
    preprocess = PreProcess()

    # Act
    result = preprocess.get_variables(config, out_dir_final)

    # Assert
    assert result.var1 == expected_var1


def test_get_userinput():
    # Arrange
    preprocess = PreProcess()
    preprocess.some_attribute = "some_value"

    # Act
    user_input = preprocess.get_userinput()

    # Assert
    assert isinstance(user_input, UserInput)
    assert user_input.some_attribute == "some_value"
