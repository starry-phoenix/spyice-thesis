class GeometrySettings:
    """Set up model geometry with two parameters, an integer geom and a float dz."""

    def __init__(self, geom: int, dz: float) -> None:
        """
        Args:
          geom (int): 1 (is a test case scenario) or  2 (is according to W3 in Buffo et al. 2018)

          dz (float): The parameter `dz` appears to be a float type. It is likely used to represent a specific value related to geometry or spatial calculations. If you need further assistance or have any specific questions about how to use this parameter in your code, feel free to ask!
        """
        self.geom = geom
        self.dz = dz
        if geom not in {1, 2}:
            error_message = "Model geometry not available"
            raise ValueError(error_message)

        self.Z = 1
        self.nc = int(self.Z / self.dz)
        self.nz = int(self.nc + 1)
