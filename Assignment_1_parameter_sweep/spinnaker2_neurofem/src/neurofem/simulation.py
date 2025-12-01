from .calc import calc_mse, calc_rmse, calc_relative_residual
from .config import NeurofemConfig, SpinnakerConfig
import numpy as np
import scipy.sparse.linalg as spla
from .problems import float_to_signed_sparse

class NeurofemSimulation:
    """
    Class representing a NeuroFEM run on SpiNNaker hardware.
    """
    def __init__(self, spinn_config: SpinnakerConfig | dict, config: NeurofemConfig | dict, matrix: np.ndarray, rhs: np.ndarray):
        """
        Create a NeuroFEM simulation instance.

        :param spinn_config: The SpiNNaker hardware configuration parameters either as a spinnaker_config object or a dictionary.
        :param config: The NeuroFEM configuration parameters as a NeuroFEMConfig object or a dictionary.
        :param matrix: The system matrix A.
        :param rhs: The right-hand side vector b.
        """
        if isinstance(config, dict):
            config = NeurofemConfig(**config)
        self.config = config
        if isinstance(spinn_config, dict):
            spinn_config = SpinnakerConfig(**spinn_config)
        self.spinn_config = spinn_config

        self.matrix = matrix
        self.rhs = rhs

        # Quantize matrix to signed integers
        self.matrix_quant, self.matrix_scale = float_to_signed_sparse(self.matrix)
        self.nmesh = len(self.rhs)

        max_cores = self.spinn_config.max_cores
        npm = self.config.npm
        self.meshes_per_core = self.nmesh // max_cores if self.nmesh % max_cores == 0 else (self.nmesh // max_cores) + 1
        self.neurons_per_core = self.meshes_per_core * npm
        self.nb_neurons = self.nmesh * npm
        self.nb_cores = self.nb_neurons // self.neurons_per_core if self.nb_neurons % self.neurons_per_core == 0 else (self.nb_neurons // self.neurons_per_core) + 1

        MAX_NEURONS_PER_CORE = 2048
        if self.neurons_per_core > MAX_NEURONS_PER_CORE:
            raise ValueError(f"Neurons per core ({self.neurons_per_core}) exceed maximum allowed ({MAX_NEURONS_PER_CORE}).")

    def _build_conns(self) -> list[list[int]]:
        """
        Generate the connection list for the SNN projection from the quantized matrix.
        :return: A list of connections in the format [pre_neuron_id, post_neuron_id, weight, delay].
        """
        row_indices, col_indices = self.matrix_quant.nonzero()
        nnz = len(row_indices)
        conns = np.zeros((nnz, 4), dtype=int)

        for idx, (i_mesh, j_mesh) in enumerate(zip(row_indices, col_indices)):
            value = self.matrix_quant[i_mesh, j_mesh]
            conns[idx] = [
                j_mesh * self.config.npm,
                i_mesh * self.config.npm,
                value,
                0
            ]
        return conns.tolist()

    def _build_neuron_params(self) -> dict:
        """
        Generate the neuron parameters dictionary for the SNN population.
        :return: A dictionary of neuron parameters.
        """
        neuron_params = {
            "gb": [self.rhs[i] * self.config.gamma for i in range(self.nmesh)],
            "threshold": self.config.theta,
            "scale": self.matrix_scale * (self.config.gamma ** 2),
            "dt": self.config.dt,
            "gamma": self.config.gamma,
            "lambda_d": self.config.lambda_d,
            "lambda_v": self.config.lambda_v,
            "k_p": self.config.k_p,
            "k_i": self.config.k_i,
            "sigma": self.config.sigma,
            "steady_state": self.config.steady_state,
        }
        return neuron_params

    def _get_solution(self, x_means, timesteps: int) -> np.ndarray:
        """
        ToDo: document

        :param timesteps: The number of timesteps the simulation was run for.
        :return The solution vector reconstructed from the neuron mean firing rates.
        """
        solution = []
        for i in range(self.nmesh):
            index = i % self.meshes_per_core + (i // self.meshes_per_core) * self.config.npm * self.meshes_per_core
            r = x_means[index] / (timesteps * self.config.steady_state + 1)
            solution.append(r)
        return np.array(solution)

    def run(self, timesteps: int = 5e+4, sys_tick_in_s: float = 1e-3, mapping_only = False) -> np.ndarray:
        """
        Run the NeuroFEM simulation on SpiNNaker hardware and retrieve the solution.
        The sys_tick_in_s parameter should not be set too low, as then too many neuron spikes may be lost.

        :param timesteps: Number of simulation timesteps to run.
        :param sys_tick_in_s: System tick duration in seconds.
        :param mapping_only: If True, only perform the mapping without running.
        :return: The solution vector as a NumPy array.
        :rtype: np.ndarray
        """
        from spinnaker2 import snn, hardware
        neuron_params = self._build_neuron_params()
        pop = snn.Population(size=self.nb_neurons, neuron_model="neurofem_2048", params=neuron_params, name=f"pop-neurofem-{id(self)}", record=["x_mean"])
        pop.set_max_atoms_per_core(self.neurons_per_core)

        conns = self._build_conns()
        proj = snn.Projection(pre=pop, post=pop, connections=conns, name=f"proj-neurofem-{id(self)}")

        net = snn.Network(f"net-neurofem-{id(self)}")
        net.add(pop, proj)

        # Run the network on SpiNNaker hardware
        hw = self.spinn_config.get_hardware()
        hw.run(net, timesteps, debug=False, sys_tick_in_s=sys_tick_in_s, mapping_only=mapping_only)

        x_mean = pop.get_x_mean()
        solution = self._get_solution(x_mean, timesteps)

        return solution
