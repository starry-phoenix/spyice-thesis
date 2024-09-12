from __future__ import annotations
import warnings

import numpy as np
from scipy.sparse import linalg

from src.spyice.parameters.user_input import UserInput as ui
from src.spyice.coefficients import update_coefficients
from src.spyice.rhs import apply_boundary_condition, correct_for_brine_movement

top_temp = ui.temperature_top_type
# TODO: Correct analytical solution
"""Solves for X in Advection- Diffusion- Equation of the form : a * (dU/dt) +  b * (dU/dz) + d/dz(c * dU/dz) + d * (dW/dt) = 0 """


class AdvectionDiffusion:
    """Class representing an advection-diffusion model."""

    def __init__(
        self,
        argument,
        X,
        source,
        X_initial,
        W,
        W_initial,
        w,
        dt,
        dz,
        nz,
        t_passed,
        S_IC,
        Stefan=False,
        Buffo=False,
        bc_neumann=None,
    ):
        """
        Args:
            argument (str): The argument for the advection-diffusion equation. Must be either "temperature" or "salinity".
            X (float): The value of X.
            source (float): The source value.
            X_initial (float): The initial value of X.
            W (float): The value of W.
            W_initial (float): The initial value of W.
            w (float): The value of w.
            dt (float): The time step size.
            dz (float): The spatial step size.
            nz (int): The number of spatial steps.
            t_passed (float): The time passed.
            S_IC (float): The initial condition for S.
            Stefan (bool, optional): Whether to use Stefan condition. Defaults to False.
            Buffo (bool, optional): Whether to use Buffo condition. Defaults to False.
            bc_neumann (float, optional): The Neumann boundary condition. Defaults to None.

        Raises:
            AssertionError: If the argument is not "temperature" or "salinity".
            AssertionError: If both Stefan and Buffo conditions are set to True.
        """

        assert argument in ["temperature", "salinity"], "invalid input for argument"
        self.argument = argument
        self.source = source
        self.Stefan = Stefan
        self.Buffo = Buffo
        self.temp_grad = bc_neumann
        assert (
            self.Stefan is not True or self.Buffo is not True
        ), "Stefan and Buffo cannot be True at the same time"
        self.nz = nz
        self.dt = dt
        self.dz = dz
        self.X_initial = X_initial
        self.w = w
        self.t_passed = t_passed
        self.S_IC = S_IC
        self.S_bc_top = "Dirichlet"
        self.top_temp = top_temp  # OR "T_const_250" OR "T_const_260" OR "T_const_265" OR "T_W3" OR "Stefan"
        ### Compute a, b, c, d
        self.Delta_W = np.zeros(nz)
        [self.a, self.b, self.c, self.d] = update_coefficients(
            argument, X_initial, w, W, nz, S_IC
        )  #
        self.Delta_W = W - W_initial  # neg: freezing pos: melting and W is phi_k
        ### Initialize
        self.main_A = np.zeros(nz, dtype=np.float64)
        self.lower_A = np.zeros(nz - 1, dtype=np.float64)
        self.upper_A = np.zeros(nz - 1, dtype=np.float64)
        self.A = np.zeros([nz, nz], dtype=np.float64)
        self.X_new = np.zeros(nz, dtype=np.float64)
        self.factor1 = self.factor_1(argument, self.a, self.c, dt, dz, nz)
        if argument == "salinity":
            self.factor1_plus, self.factor1_minus = self.factor1
        self.factor2 = self.factor_2(self.a, self.b, dt, dz, nz)
        self.factor3 = self.factor_3(self.a, self.d, nz)

    def TDMAsolver(self, a, b, c, d):
        """Solves a tridiagonal matrix system using the TDMA (Thomas algorithm).
        Args:
            a (array-like): Coefficients of the sub-diagonal elements.
            b (array-like): Coefficients of the diagonal elements.
            c (array-like): Coefficients of the super-diagonal elements.
            d (array-like): Right-hand side vector.
        Returns:
            array-like: Solution vector.
        References:
            - Tridiagonal matrix algorithm: http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
            - TDMA (Thomas algorithm): http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
        """
        # TDMA solver, a b c d can be NumPy array type or Python list type.
        # refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        # and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)

        nf = len(d)  # number of equations
        ac, bc, cc, dc = map(np.array, (a, b, c, d))  # copy arrays
        for it in range(1, nf):
            mc = ac[it - 1] / bc[it - 1]
            bc[it] = bc[it] - mc * cc[it - 1]
            dc[it] = dc[it] - mc * dc[it - 1]

        xc = bc
        xc[-1] = dc[-1] / bc[-1]

        for il in range(nf - 2, -1, -1):
            xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

        return xc

    # ---------WORKING CODE BUFFO MATLAB-----------------
    def Buffosolver(self, a, b, c, f):
        """Solves a tridiagonal linear system using the Thomas algorithm.
        Args:
            a (np.ndarray): Lower diagonal of the tridiagonal matrix.
            b (np.ndarray): Upper diagonal of the tridiagonal matrix.
            c (np.ndarray): Main diagonal of the tridiagonal matrix.
            f (np.ndarray): Right-hand side vector.
        Returns:
            np.ndarray: Solution vector of the linear system.
        """

        a, b, c, f = map(
            np.array, (a, b, c, f)
        )  # lower, upper, main diagonal, RHS-vector
        assert (
            np.round(np.abs(a) + np.abs(b), 2) < np.round(np.abs(c[:-1]))
        ).all(), (
            f"Matrix is not diagonally dominant, TDMA fails at {self.t_passed}s"
        )  # checks if c > a + b upto 2nd order

        n = len(f)
        alpha = np.zeros(n)
        beta = np.zeros(n)
        x = np.zeros(n)
        a = np.append(0, a)
        b = np.append(b, 0)
        alpha[0] = b[0] / c[0]
        beta[0] = f[0] / c[0]

        for i in range(1, n):
            alpha[i] = b[i] / (-a[i] * alpha[i - 1] + c[i])
            beta[i] = (f[i] - a[i] * beta[i - 1]) / (-a[i] * alpha[i - 1] + c[i])

        x[n - 1] = beta[n - 1]

        for i in reversed(range(n - 1)):
            x[i] = -alpha[i] * x[i + 1] + beta[i]

        return x

    def set_up_tridiagonal(self):
        """Set up the tridiagonal matrix for Salinity and Temperature.
        This method sets up the tridiagonal matrix for Salinity and Temperature calculations. It initializes the main diagonal (main_A), upper diagonal (upper_A), and lower diagonal (lower_A) of the matrix based on the given arguments.
        Args:
            self: The object instance.
        Returns:
            None
        """

        ### Set up tridiagonal matrix LHS for Salinity and Temperature
        for i in range(self.nz):
            if self.argument == "salinity":
                self.main_A[i] = 1 + self.factor1_plus[i] + self.factor1_minus[i]
                if i < self.nz - 1:
                    self.upper_A[i] = -self.factor1_minus[i]
                if i > 0:
                    self.lower_A[i] = -self.factor1_plus[i]
            else:
                self.main_A[i] = 2 * self.factor1[i] + 1.0
                if i < self.nz - 1:
                    self.upper_A[i] = -self.factor1[i]
                if i > 0:
                    self.lower_A[i - 1] = -self.factor1[i]

        if self.argument == "temperature":
            if self.Stefan is True:
                self.lower_A[-1] = 0.0
                self.main_A[-1] = 1.0
                if self.temp_grad is not None:
                    self.upper_A[0] = self.upper_A[0] - self.factor1[0]
                    self.main_A[0] = self.main_A[0]
                else:
                    self.main_A[0] = 1.0
                    self.upper_A[0] = 0.0

            elif self.Buffo is True:
                self.upper_A[0] = self.upper_A[0] - self.factor1[0]
            else:
                # non-pragmatic dirichlet RB like Buffo
                self.upper_A[0] = -1 * self.factor1[0]
                self.main_A[0] = 2 * self.factor1[0] + 1
                self.main_A[-1] = 2 * self.factor1[-1] + 1
                self.lower_A[-1] = -1 * self.factor1[-1]

        # self.main_A[0] = self.main_A[0]  #
        # self.main_A[-1] = 0
        ### Set up tridiagonal Matrix A and solve for new u
        self.A = np.zeros([self.nz, self.nz])
        self.A = (
            np.diag(self.main_A, k=0)
            + np.diag(self.lower_A, k=-1)
            + np.diag(self.upper_A, k=1)
        )

        # return self.lower_A, self.main_A, self.upper_A

    def unknowns_matrix(self):
        """Calculates the unknowns matrix for the advection-diffusion model.
        Returns:
            tuple: A tuple containing the following elements:
                - X_new (ndarray): The solution vector for the unknowns.
                - X_wind (ndarray): The solution vector corrected for brine movement.
                - dt (float): The time step used in the calculation.
        """
        # Rest of the code...

        self.set_up_tridiagonal()

        B = apply_boundary_condition(
            self.argument,
            self.X_initial,
            self.source,
            self.factor1,
            self.factor3,
            self.a,
            self.Delta_W,
            self.w,
            self.nz,
            self.t_passed,
            self.S_IC,
            self.top_temp,
            self.Stefan,
            self.Buffo,
            self.temp_grad,
        )

        X_wind = correct_for_brine_movement(
            self.argument,
            self.X_initial,
            self.w,
            self.t_passed,
            self.nz,
            self.S_IC,
            self.top_temp,
        )

        if self.Buffo is True:
            X_new = self.Buffosolver(
                self.lower_A, self.upper_A, self.main_A, B
            )  # input: lower, upper, main diagonal, RHS-vector
        else:
            warnings.filterwarnings("ignore")
            X_new = linalg.spsolve(self.A, B)
        return X_new, X_wind

    def factor_1(self, argument, a, c, dt, dz, nz):
        """Factor 1 and avoid zero divison error

        Args:
            argument (str): The argument value.
            a (numpy.ndarray): The array containing values for 'a'.
            c (numpy.ndarray): The array containing values for 'c'.
            dt (float): The value of 'dt'.
            dz (float): The value of 'dz'.
            nz (int): The value of 'nz'.
        Returns:
            numpy.ndarray or list: The calculated factor(s) based on the given argument.
                If the argument is not "salinity", returns a numpy.ndarray.
                If the argument is "salinity", returns a list containing two numpy.ndarrays.

        """
        const1 = dt / (dz**2)
        self.factor1 = np.zeros(nz, dtype=np.float64)
        self.factor1[np.nonzero(a)] = const1 * (c[np.nonzero(a)] / a[np.nonzero(a)])

        if argument != "salinity":
            return self.factor1
        factor1_plus = (self.factor1 + np.roll(self.factor1, 1)) / 2
        factor1_plus[0] = self.factor1[0]
        factor1_minus = (self.factor1 + np.roll(self.factor1, -1)) / 2
        factor1_minus[-1] = self.factor1[-1]
        return [factor1_plus, factor1_minus]

    def factor_2(self, a, b, dt, dz, nz):
        """Calculates the factor2 array for advection-diffusion model.
        Args:
            a (numpy.ndarray): Array containing values for 'a'.
            b (numpy.ndarray): Array containing values for 'b'.
            dt (float): Time step.
            dz (float): Spatial step.
            nz (int): Number of elements in the array.
        Returns:
            numpy.ndarray: Array containing the calculated factor2 values.
        """

        factor2 = np.zeros(nz, dtype=np.float64)
        const2 = dt / dz
        factor2[np.nonzero(a)] = const2 * (b[np.nonzero(a)] / a[np.nonzero(a)])
        return factor2

    def factor_3(self, a, d, nz):
        """Factor 3 and avoid zero division error.

        Args:
            a (numpy.ndarray): Array containing values.
            d (numpy.ndarray): Array containing values.
            nz (int): Number of elements.
        Returns:
            numpy.ndarray: Array containing factor 3 values.

        """
        factor3 = np.zeros(nz, dtype=np.float64)
        factor3[np.nonzero(a)] = d[np.nonzero(a)] / a[np.nonzero(a)]
        return factor3


# %%
