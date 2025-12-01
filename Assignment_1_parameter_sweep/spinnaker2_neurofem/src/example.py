from argparse import ArgumentParser
from neurofem import generate_poisson_unitdisk_variable_f, float_to_signed_sparse, NeurofemSimulation, NeurofemConfig, SpinnakerConfig
from neurofem.calc import calc_mse, calc_rmse, calc_relative_residual

def main():
    parser = ArgumentParser(description="Run NeuroFEM simulation on SpiNNaker hardware.")
    parser.add_argument("-g", type=float, default=50, help="Gamma parameter for NeuroFEM configuration.")
    parser.add_argument("-dt", type=float, default=None, help="Time step (dt) for NeuroFEM configuration.")
    parser.add_argument("-t", type=int, default=100000, help="Number of timesteps")
    parser.add_argument("-st", type=float, default=1e-3, help="Simulation time step (dt) for SpiNNaker hardware -- (sys_tick_in_s).")
    args = parser.parse_args()

    gamma = args.g
    dt = args.dt
    timesteps = args.t
    sys_tick_in_s = args.st

    matrix, rhs, _ = generate_poisson_unitdisk_variable_f()
    spinn_config = SpinnakerConfig(s2_ip="192.168.0.6")
    neurofem_config = NeurofemConfig(gamma=gamma, dt=dt)
    simulation = NeurofemSimulation(spinn_config, neurofem_config, matrix, rhs)

    solution = simulation.run(timesteps=timesteps, sys_tick_in_s=sys_tick_in_s)
    mse = calc_mse(matrix, rhs, solution)
    rmse = calc_rmse(matrix, rhs, solution)
    rel_residual = calc_relative_residual(matrix, rhs, solution)

    print(f"MSE: {mse}")
    print(f"RMSE: {rmse}")
    print(f"Relative Residual: {rel_residual}")

if __name__ == "__main__":
    main()
