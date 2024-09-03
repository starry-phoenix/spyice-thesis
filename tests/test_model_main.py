import pytest
from spyice.main_process import MainProcess


@pytest.fixture
def main_process():
    config = {
        "consts": {"constants": "real"},
        "dt": {"dt": 47.0},
        "S_IC": {"S_IC": "S34"},
        "iter_max": {"iter_max": 1000},
    }
    project_path = "~/RWTHMSc/Sem2/MBD_Parttime/Sea-Ice-Model-Sneha/spyice/spyice/tests"
    hyd_output_dir = "./outputs/real-S34-iter_max1000/0"

    return MainProcess(config, hyd_output_dir=hyd_output_dir, project_path=project_path)


def test_main_process_creation(main_process):
    # Arrange
    # config = {
    #     "consts": {"constants": "real"},
    #     "dt": {"dt": 47.0},
    #     "S_IC": {"S_IC": "S34"},
    #     "iter_max": {"iter_max": 1000},
    # }
    # project_path = "~/RWTHMSc/Sem2/MBD_Parttime/Sea-Ice-Model-Sneha/spyice/spyice/tests"
    # hyd_output_dir = "./outputs/real-S34-iter_max1000/0"

    # # Act
    # main_process = MainProcess(
    #     config, hyd_output_dir=hyd_output_dir, project_path=project_path
    # )

    # Assert

    assert hasattr(main_process, "project_path")
    assert hasattr(main_process, "config")
    assert hasattr(main_process, "run_model")
    assert main_process.config == {
        "consts": {"constants": "real"},
        "dt": {"dt": 47.0},
        "S_IC": {"S_IC": "S34"},
        "iter_max": {"iter_max": 1000},
    }
    assert (
        main_process.project_path
        == "~/RWTHMSc/Sem2/MBD_Parttime/Sea-Ice-Model-Sneha/spyice/spyice/tests"
    )
