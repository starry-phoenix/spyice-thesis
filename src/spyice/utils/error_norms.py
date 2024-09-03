import numpy as np


class ErrorNorms:
    """Defines functions for calculating error norms between numerical and analytical values."""

    def __init__(self, numerical_values, analytical_values) -> None:
        """
        Args:
            numerical_values (list): A list of numerical values.
            analytical_values (list): A list of analytical values.

        Attributes:
            numerical_values (list): A list of numerical values.
            analytical_values (list): A list of analytical values.
            iteration_count (int): The count of iterations.

        """
        self.numerical_values = numerical_values
        self.analytical_values = analytical_values
        self.iteration_count = 0

    def numerical_analytical_diff(self):
        """Calculates the absolute difference between the analytical values and the numerical values.

        Args:
            None

        Returns:
            numpy.ndarray: The array containing the absolute differences between the analytical values and the numerical values.
        """

        if len(self.analytical_values.shape) == 1:
            self.iteration_count, self.array_len = 1, self.analytical_values.shape[0]
            self.numerical_analytical_diff = abs(
                self.analytical_values - self.numerical_values
            )
        else:
            [self.iteration_count, self.array_len] = self.analytical_values.shape
            self.numerical_analytical_diff = np.array(
                [
                    abs(self.analytical_values[i] - self.numerical_values[i])
                    for i in range(self.iteration_count)
                ]
            )
        return self.numerical_analytical_diff

    def one_norm(self, numerical_analytical_diff):
        """Calculates the one norm of the numerical-analytical difference.

        Args:
            numerical_analytical_diff (ndarray): The numerical-analytical difference.
        Returns:
            ndarray or float: The one norm of the numerical-analytical difference. If the analytical values are 1D, a float is returned. Otherwise, an ndarray is returned with the one norm for each iteration.
        Raises:
            None
        """

        return (
            np.sum(numerical_analytical_diff) / self.array_len
            if len(self.analytical_values.shape) == 1
            else np.array(
                [
                    np.sum(numerical_analytical_diff[i]) / self.array_len
                    for i in range(self.iteration_count)
                ]
            )
        )

    def two_norm(self, numerical_analytical_diff):
        """Calculates the two-norm error of the numerical-analytical difference.

        Args:
            numerical_analytical_diff (ndarray): The numerical-analytical difference.
        Returns:
            ndarray or float: The two-norm error. If the analytical values are 1D, a float is returned.
            If the analytical values are 2D, an ndarray is returned with the two-norm error for each iteration.
        """

        return (
            np.sqrt(np.sum(numerical_analytical_diff**2)) / self.array_len
            if len(self.analytical_values.shape) == 1
            else np.array(
                [
                    np.sqrt(np.sum(numerical_analytical_diff[i] ** 2)) / self.array_len
                    for i in range(self.iteration_count)
                ]
            )
        )

    def infinity_norm(self, numerical_analytical_diff):
        """Calculates the infinity norm of the given numerical-analytical difference.

        Args:
            numerical_analytical_diff (array-like): The numerical-analytical difference.
        Returns:
            float or numpy.ndarray: The infinity norm of the difference. If the analytical values
            are 1-dimensional, a single float value is returned. Otherwise, a numpy array
            containing the infinity norm for each iteration is returned.
        """

        return (
            max(numerical_analytical_diff)
            if len(self.analytical_values.shape) == 1
            else np.array(
                [max(numerical_analytical_diff[i]) for i in range(self.iteration_count)]
            )
        )
