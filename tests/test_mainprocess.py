import pytest
from unittest.mock import patch, MagicMock

# Import the function to be tested
from spyice.main_process import MainProcess


@pytest.mark.parametrize(
    "config, out_dir_final, preprocess_vars, userinput, results, analysis, expected_plot_call",
    [
        # Happy path test case
        (
            {"param1": "value1"},
            "output/dir",
            {"var1": "data1"},
            {"input1": "data2"},
            MagicMock(t_k_diff="diff1", t_stefan_diff="diff2"),
            {"error1": "result1"},
            (
                {"input1": "data2"},
                MagicMock(t_k_diff="diff1", t_stefan_diff="diff2"),
                {"error1": "result1"},
            ),
        ),
        # Edge case with empty config and out_dir_final
        (
            {},
            "",
            {},
            {},
            MagicMock(t_k_diff="", t_stefan_diff=""),
            {},
            ({}, MagicMock(t_k_diff="", t_stefan_diff=""), {}),
        ),
        # Error case with None values
        (None, None, None, None, None, None, (None, None, None)),
    ],
    ids=["happy_path", "edge_case_empty_values", "error_case_none_values"],
)
def test_run_model(
    config,
    out_dir_final,
    preprocess_vars,
    userinput,
    results,
    analysis,
    expected_plot_call,
):
    # Arrange
    self = MagicMock()
    self.config = config
    self.out_dir_final = out_dir_final

    with patch(
        "spyice.preprocess.model_preprocess.PreProcess.get_variables",
        return_value=preprocess_vars,
    ) as mock_get_variables, patch(
        "spyice.preprocess.model_preprocess.PreProcess.get_userinput",
        return_value=userinput,
    ) as mock_get_userinput, patch(
        "spyice.model.src.run_process.SeaIceModelClass.get_results",
        return_value=results,
    ) as mock_get_results, patch(
        "spyice.postprocess.Analysis.get_error_results", return_value=analysis
    ) as mock_get_error_results, patch.object(self, "plot_model") as mock_plot_model:
        # Act
        MainProcess.run_model(self)

        # Assert
        mock_get_variables.assert_called_once_with(config, out_dir_final)
        mock_get_userinput.assert_called_once()
        mock_get_results.assert_called_once_with(preprocess_vars)
        mock_get_error_results.assert_called_once_with(
            t_k_diff=results.t_k_diff, t_stefan_diff=results.t_stefan_diff
        )
        mock_plot_model.assert_called_once_with(*expected_plot_call)
