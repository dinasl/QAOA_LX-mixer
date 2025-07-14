from multiprocessing import Pool

#import tikzplotlib

import itertools
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving plots
from matplotlib import pyplot as plt
import numpy as np
from random import sample
import time
from Mixer import *

# # Importing necessary classes from Mixer module
# import sys
# sys.path.append(r"C:\Users\sanne\LogicalXMixer")
# from Mixer_Franz import MixerFranz

# m = Mixer(B=["100", "010", "001", "110", "101", "011", "111"])

# print(m)
class Worker:

    """n = dim(hilbert space)"""

    def __init__(self, n, checkmixer=False):
        self.n = n  # Store n for later use
        self.all_states = ["".join(i) for i in itertools.product("01", repeat=n)]
        self.checkmixer = checkmixer

    def sample_B(self, m):
        """m = |B|"""
        binary_strings = sample(self.all_states, m)
        # Convert binary strings to integers for LXMixer
        return [int(binary_str, 2) for binary_str in binary_strings]

    def get_costs(self, B):
        nL = self.n  # Use the stored number of qubits
        
        # Finding the cnot cost for a given B and nL
        # Note: keeping loop structure for future chain mixer comparison
        for chain in [True, False]:
            if chain:
                # Future: chain mixer implementation
                continue
            else:
                # The LXMixer class
                mixer = LXMixer(B, nL)
                mixer.compute_family_of_valid_graphs()
                mixer.compute_all_orbits()
                stabilizer = Stabilizer(B=mixer.B, n=mixer.nL, orbit_dictionary=mixer.orbits)
                stabilizer.compute_minimal_generating_sets()
                stabilizer.compute_projector_stabilizers()
                mixer.compute_costs()

                # Find the best mixer configuration
                best_Xs, best_Zs, best_cost = mixer.find_best_mixer()
                cnots_optimal = best_cost
                break  # Exit loop after processing LXMixer

        print("cnots_optimal", cnots_optimal)
        return cnots_optimal#cnots, cnots_chain, num_optimal_solutions


def saveResult(res):
    global cnots_optimal #cnots, cnots_reduced, cnots_chain, cnots_chain_reduced, num_optimal_solutions
    cnots_optimal.append(res)
    # cnots.append(res[0])
    # cnots_reduced.append(res[1])
    # cnots_chain.append(res[2])
    # cnots_chain_reduced.append(res[3])
    # num_optimal_solutions.append(res[4])


def plot(m_list, mean_cnots, var_cnots, min_cnots, max_cnots, color, col, style, text):
    plt.plot(m_list, mean_cnots, style+"-"+col, label=r"mean "+text)
    plt.fill_between(
        m_list,
        mean_cnots - np.sqrt(var_cnots),
        mean_cnots + np.sqrt(var_cnots),
        color=color,
        alpha=0.35,
        label=r"$1\sigma$",
    )
    plt.plot(m_list, min_cnots, style+":k")
    plt.plot(m_list, max_cnots, style+":k", label=r"min/max "+text)
    # plt.plot(m_list, mean_cnots_opt-np.sqrt(var_cnots_opt), '-g', alpha=.35)
    # plt.plot(m_list, mean_cnots_opt+np.sqrt(var_cnots_opt), '-g', alpha=.35,label="std dev")
    plt.fill_between(m_list, min_cnots, max_cnots, color="black", alpha=0.15)
    # plt.plot(m_list, time_list / max(time_list), label = f"n = {n}", color = "black")
    plt.legend()
    plt.xlim(min(m_list), max(m_list))
    plt.xlabel(r"$|B|$")
    # plt.ylabel(r"$t / t_{max}$")
    plt.ylabel("#CNOTS")
    plt.grid()


def main(n, num_samples=100):
    print("n=", n)
    worker = Worker(n, checkmixer=False)

    # TODO this is a new loop for debugging... Test the worker directly first
    print("Testing worker.get_costs directly...")
    try:
        test_B = worker.sample_B(2)
        print(f"Test B: {test_B}")
        test_result = worker.get_costs(test_B)
        print(f"Test result: {test_result}")
    except Exception as e:
        print(f"Error in direct test: {e}")
        import traceback
        traceback.print_exc()
        return

    # m = |B|
    m_list = list(range(2, 2**n + 1))

    # run
    # time_list = np.zeros(len(m_list))

    min_cnots_optimal = np.zeros(len(m_list))
    max_cnots_optimal = np.zeros(len(m_list))
    mean_cnots_optimal = np.zeros(len(m_list))
    var_cnots_optimal = np.zeros(len(m_list))

    # min_cnots = np.zeros(len(m_list))
    # max_cnots = np.zeros(len(m_list))
    # mean_cnots = np.zeros(len(m_list))
    # var_cnots = np.zeros(len(m_list))

    # min_cnots_reduced = np.zeros(len(m_list))
    # max_cnots_reduced = np.zeros(len(m_list))
    # mean_cnots_reduced = np.zeros(len(m_list))
    # var_cnots_reduced = np.zeros(len(m_list))

    # min_cnots_chain = np.zeros(len(m_list))
    # max_cnots_chain = np.zeros(len(m_list))
    # mean_cnots_chain = np.zeros(len(m_list))
    # var_cnots_chain = np.zeros(len(m_list))

    # min_cnots_chain_reduced = np.zeros(len(m_list))
    # max_cnots_chain_reduced = np.zeros(len(m_list))
    # mean_cnots_chain_reduced = np.zeros(len(m_list))
    # var_cnots_chain_reduced = np.zeros(len(m_list))

    # min_num_optimal_solutions = np.zeros(len(m_list))
    # max_num_optimal_solutions = np.zeros(len(m_list))
    # mean_num_optimal_solutions = np.zeros(len(m_list))
    # var_num_optimal_solutions = np.zeros(len(m_list))

    for i, m in enumerate(m_list):
        global cnots_optimal #cnots, cnots_reduced, cnots_chain, cnots_chain_reduced, num_optimal_solutions
        # cnots = []
        # cnots_reduced = []
        # cnots_chain = []
        # cnots_chain_reduced = []
        # num_optimal_solutions = []
        cnots_optimal = []
        pool = Pool()
        results = []
        for j in range(num_samples):
            result = pool.apply_async(
                worker.get_costs, args=(worker.sample_B(m),), callback=saveResult
            )
            # try:
            #     deb.get()
            # except Exception as e:
            #     print("Exception in worker.run:", e)
            #     traceback.print_exc()
            results.append(result)
        pool.close()
        pool.join()
        
        # TODO this is a new debugging loop (until min_cnots_optimal)... Check for exceptions in the results
        for j, result in enumerate(results):
            try:
                result.get()  # This will raise any exception that occurred
            except Exception as e:
                print(f"Exception in worker {j}: {e}")
                import traceback
                traceback.print_exc()
        
        print(f"Collected {len(cnots_optimal)} results for m={m}")
        if len(cnots_optimal) == 0:
            print("No results collected - all tasks failed!")
            continue  # Skip this m value

        min_cnots_optimal[i] = np.min(cnots_optimal)
        max_cnots_optimal[i] = np.max(cnots_optimal)
        mean_cnots_optimal[i] = np.mean(cnots_optimal)
        var_cnots_optimal[i] = np.var(cnots_optimal)

        # min_cnots[i] = np.min(cnots)
        # max_cnots[i] = np.max(cnots)
        # mean_cnots[i] = np.mean(cnots)
        # var_cnots[i] = np.var(cnots)

        # min_cnots_reduced[i] = np.min(cnots_reduced)
        # max_cnots_reduced[i] = np.max(cnots_reduced)
        # mean_cnots_reduced[i] = np.mean(cnots_reduced)
        # var_cnots_reduced[i] = np.var(cnots_reduced)

        # min_cnots_chain[i] = np.min(cnots_chain)
        # max_cnots_chain[i] = np.max(cnots_chain)
        # mean_cnots_chain[i] = np.mean(cnots_chain)
        # var_cnots_chain[i] = np.var(cnots_chain)

        # min_cnots_chain_reduced[i] = np.min(cnots_chain_reduced)
        # max_cnots_chain_reduced[i] = np.max(cnots_chain_reduced)
        # mean_cnots_chain_reduced[i] = np.mean(cnots_chain_reduced)
        # var_cnots_chain_reduced[i] = np.var(cnots_chain_reduced)

        # min_num_optimal_solutions[i] = np.min(num_optimal_solutions)
        # max_num_optimal_solutions[i] = np.max(num_optimal_solutions)
        # mean_num_optimal_solutions[i] = np.mean(num_optimal_solutions)
        # var_num_optimal_solutions[i] = np.var(num_optimal_solutions)

        print(int(100 * (i + 1) / len(m_list)), "%")

    fig = plt.figure()

    # plot(
    #     m_list,
    #     mean_cnots_reduced,
    #     var_cnots_reduced,
    #     min_cnots_reduced,
    #     max_cnots_reduced,
    #     color="green",
    #     col="g",
    #     style="x",
    #     text=", optimal reduced"
    # )
    # #tikzplotlib.save("statistics_cnots_n" + str(n) + "_1.tex")
    # plt.savefig("stat_n" + str(n) + "_1.png")

    # plot(
    #     m_list,
    #     mean_cnots,
    #     var_cnots,
    #     min_cnots,
    #     max_cnots,
    #     color="blue",
    #     col="b",
    #     style="o",
    #     text=", optimal"
    # )
    # #tikzplotlib.save("statistics_cnots_n" + str(n) + "_2.tex")
    # plt.savefig("stat_n" + str(n) + "_2.png")

    # plot(
    #     m_list,
    #     mean_cnots_chain,
    #     var_cnots_chain,
    #     min_cnots_chain,
    #     max_cnots_chain,
    #     color="red",
    #     col="r",
    #     style="+",
    #     text=", $T_{\Delta}$"
    # )
    # #tikzplotlib.save("statistics_cnots_n" + str(n) + "_3.tex")
    # plt.savefig("stat_n" + str(n) + "_3.png")

    # plot(
    #     m_list,
    #     mean_cnots_chain_reduced,
    #     var_cnots_chain_reduced,
    #     min_cnots_chain_reduced,
    #     max_cnots_chain_reduced,
    #     color="yellow",
    #     col="y",
    #     style="1",
    #     text=", $T_{\Delta}$, reduced"
    # )
    # #tikzplotlib.save("statistics_cnots_n" + str(n) + "_4.tex")
    # plt.savefig("stat_n" + str(n) + "_4.png")

    # plt.clf()

    # plot(
    #     m_list,
    #     mean_num_optimal_solutions,
    #     var_num_optimal_solutions,
    #     min_num_optimal_solutions,
    #     max_num_optimal_solutions,
    #     color="black",
    #     col="k",
    #     style="1",
    #     text=", num opt sol"
    # )
    # #tikzplotlib.save("statistics_num_opt_sol_n" + str(n) + ".tex")
    # plt.savefig("stat_n" + str(n) + "_num_opt_sol.png")

    # plt.clf()

    plot(
        m_list,
        mean_cnots_optimal,
        var_cnots_optimal,
        min_cnots_optimal,
        max_cnots_optimal,
        color="green",
        col="g",
        style="x",
        text=", optimal reduced"
    )
    #tikzplotlib.save("statistics_cnots_n" + str(n) + "_1.tex")
    plt.savefig("stat_n" + str(n) + "_1.png")
    plt.clf()

    
if __name__ == "__main__":
    pass
    # n=dimension of Hilbert space
    # for n in [3]:
    #     main(n)