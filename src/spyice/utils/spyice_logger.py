from io import StringIO

# Step 2: Assuming Hydra and logging are already configured elsewhere in your application


# Step 3: Create a logging replacement for print
class SpyiceLogger(StringIO):
    """A class that replaces the print function with a logger object."""

    def __init__(self, logger, *args, **kwargs):
        """
        Args:
            logger: The logger object used for logging.
            *args: Additional positional arguments.
            **kwargs: Additional keyword arguments.
        """
        super().__init__(*args, **kwargs)
        self.logger = logger

    def write(self, message):
        """Writes the given message to the logger if it is not empty.

        Args:
            message (str): The message to be written to the logger.
        Returns:
            None
        """

        if message.rstrip():
            self.logger.info(message)

    def flush(self):
        pass
