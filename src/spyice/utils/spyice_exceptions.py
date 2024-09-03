class InvalidPhaseError(Exception):
    """Custom exception class for SeaIceModel errors."""


class ConvergenceError(Exception):
    """Exception raised when convergence is not reached."""

    def __init__(self, message="Convergence not reached."):
        self.message = message
        super().__init__(self.message)
