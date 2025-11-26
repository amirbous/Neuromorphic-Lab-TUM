import scipy.sparse.linalg as spla
import numpy as np

class SpinnakerConfig:
    """
    This class holds configuration parameters for SpiNNaker2 hardware.
    It supports both single-chip and multi-node (48-chip) boards.
    """
    CORES_PER_CHIP: int = 148
    CHIPS_PER_BOARD: int = 48
    def __init__(self, s2_ip: str | None = None, stm_ip: str | None = None, is_multi_node_board: bool = False):
        """
        Create a new SpiNNaker2 configuration.
        When using a single-chip board, provide the s2_ip. For multi-node boards, provide the stm_ip.

        :param s2_ip: IP address of the SpiNNaker2 chip (for single-chip boards)
        :param stm_ip: IP address of the STM controller
        :param is_multi_node_board: Whether using a multi-node board (48 chips) or single-chip board
        """
        if not is_multi_node_board and s2_ip is None:
            raise ValueError("s2_ip must be provided for single-chip boards.")
        if is_multi_node_board and stm_ip is None:
            raise ValueError("stm_ip must be provided for multi-node boards.")
        self._s2_ip = s2_ip
        self._stm_ip = stm_ip
        self._is_multi_node_board = is_multi_node_board

    @property
    def s2_ip(self) -> str | None:
        """
        The IP address of the SpiNNaker2 chip (for single-chip boards).
        """
        return self._s2_ip

    @s2_ip.setter
    def s2_ip(self, value: str):
        if self.is_multi_node_board:
            raise ValueError("s2_ip should not be set for multi-node boards.")
        self._s2_ip = value

    @property
    def stm_ip(self) -> str | None:
        """
        The IP address of the STM controller (for multi-node boards).
        """
        return self._stm_ip

    @stm_ip.setter
    def stm_ip(self, value: str):
        if not self.is_multi_node_board:
            raise ValueError("stm_ip should not be set for single-chip boards.")
        self._stm_ip = value

    @property
    def is_multi_node_board(self) -> bool:
        """
        Whether using a multi-node (48-chip) board.
        """
        return self._is_multi_node_board

    @property
    def max_cores(self) -> int:
        """
        Maximum number of cores available on the board.
        """
        if self.is_multi_node_board:
            return self.CORES_PER_CHIP * self.CHIPS_PER_BOARD
        else:
            return self.CORES_PER_CHIP

    def get_hardware(self) -> "hardware.SpiNNaker2Chip | hardware.SpiNNcloud48NodeBoard":
        """
        Get the py-spinnaker2 hardware object based on the configuration.+
        Note that this requires the spinnaker2 package to be installed.

        :return: SpiNNaker2 hardware object
        """
        from spinnaker2 import hardware
        from spinnaker2.experiment_backends import BackendSettings, ExperimentBackendType
        from spinnaker2.experiment_backends.backend_settings import ROUTING
        from spinn_machine.version.version_factory import SPIN2_GEN

        settings = BackendSettings()
        settings.routing_type = ROUTING.C2C
        settings.rebuild_apps = True

        kwargs = {
            "experiment_backend_type": ExperimentBackendType.SPINNMAN2,
        }
        if self.is_multi_node_board:
            kwargs["stm_ip"] = self.stm_ip
            hw = hardware.SpiNNcloud48NodeBoard(**kwargs)
        else:
            kwargs["eth_ip"] = self.s2_ip
            hw = hardware.SpiNNaker2Chip(**kwargs)
        hw.experiment_backend.board_type = SPIN2_GEN.SPIN2_48CHIP
        hw.experiment_backend.stm_ip = self.stm_ip
        return hw

class NeurofemConfig:
    """
    A class to hold configuration parameters for NeuroFEM simulations.
    This only includes parameters specific to the NeuroFEM algorithm.
    """
    def __init__(self,
        gamma: int = 50,
        dt: float = 2**-12,
        theta: float | None = None,
        lambda_d: float= 0.2,
        lambda_v: float = 0.4,
        k_p: float = 4.0,
        k_i: float = 16.0,
        sigma: float = 0.00225,
        steady_state: float = 0.4,
    ):
        """
        Create a new NeuroFEM configuration.

        | Parameter    | Value         | Description                           | Effect |
        | ------------ | ------------- | ------------------------------------- | ------ |
        | gamma        | 50            | scaling factor for the system         |        |
        | dt           | 2^-13         | simulation time step                  |        |
        | npm          | 8             | neurons per node                      |        |
        | theta        | 0.5 * gamma^2 | neuron firing threshold               |        |
        | lambda_d     | 0.2           | damping coefficient                   |        |
        | lambda_v     | 0.4           | velocity coefficient                  |        |
        | k_p          | 4.0           | proportional gain                     |        |
        | k_i          | 16.0          | integral gain                         |        |
        | sigma        | 0.00225       | noise standard deviation              |        |
        | steady_state | 0.4           | steady state firing rate factor       |        |

        :param gamma: Scaling factor for the system
        :param dt: Simulation time step ...
        :param theta: Neuron firing threshold, if None, computed as 0.5 * gamma^2
        :param lambda_d: Damping coefficient
        :param lambda_v: Velocity coefficient
        :param k_p: Proportional gain
        :param k_i: Integral gain
        :param sigma: Noise standard deviation
        :param steady_state: Steady state firing rate factor
        """
        self.gamma = gamma
        self.dt = dt
        self.theta = theta
        self.lambda_d = lambda_d
        self.lambda_v = lambda_v
        self.k_p = k_p
        self.k_i = k_i
        self.sigma = sigma
        self.steady_state = steady_state

        self.npm = 8 # Neurons per mesh node
