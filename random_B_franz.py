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

# Importing necessary classes from Mixer module
import sys
sys.path.append(r"C:\Users\sanne\LogicalXMixer")
from Mixer_Franz import MixerFranz


class Worker:

    """n = dim(hilbert space)"""

    def __init__(self, n, checkmixer=False):
        self.n = n  # Store n for later use
        self.all_states = ["".join(i) for i in itertools.product("01", repeat=n)]
        self.checkmixer = checkmixer

    def sample_B(self, m):
        """m = |B|"""
        binary_strings = sample(self.all_states, m)
        binary_integers = [int(binary_str, 2) for binary_str in binary_strings]
        # Return both formats for flexibility
        return binary_strings, binary_integers
    

    def get_costs(self, B_strings, B_integers):
        nL = self.n  # Use the stored number of qubits
        
        # Finding the cnot cost for both Franz mixer (strings) and LXMixer (integers)
        results = {}
        
        # Franz mixer (uses strings) - following the exact pattern from random_B_franz_string.py
        try:
            # Replicate the exact pattern from random_B_franz_string.py
            cnots = None
            cnots_reduced = None
            
            for chain in [True, False]:
                for reduced in [True, False]:
                    mixer_franz = MixerFranz(B_strings, reduced=reduced)
                    
                    if chain:
                        mixer_franz.get_chain_mixer()
                    else:
                        mixer_franz.get_best_mixer_commuting_graphs()
                    
                    if chain and reduced:
                        cnots_chain_reduced = mixer_franz.solution_chain_reduced_cost
                    elif chain and (not reduced):
                        cnots_chain = mixer_franz.solution_chain_cost
                    elif (not chain) and reduced:
                        cnots_reduced = mixer_franz.solution_reduced_cost
                    else:
                        cnots = mixer_franz.solution_cost
            
            # Use the "optimal" case (not chain, not reduced)
            results['franz'] = cnots
            print(f"Franz mixer - optimal cost: {cnots}")
            
        except Exception as e:
            print(f"Error with Franz mixer: {e}")
            import traceback
            traceback.print_exc()
            results['franz'] = float('inf')
        
        # LXMixer (uses integers)
        try:
            mixer = LXMixer(B_integers, nL)
            mixer.compute_family_of_valid_graphs()
            mixer.compute_all_orbits()
            stabilizer = Stabilizer(B=mixer.B, n=mixer.nL, orbit_dictionary=mixer.orbits)
            stabilizer.compute_minimal_generating_sets()
            stabilizer.compute_projector_stabilizers()
            mixer.compute_costs()
            
            best_Xs, best_Zs, best_cost = mixer.find_best_mixer()
            results['lxmixer'] = best_cost
        except Exception as e:
            print(f"Error with LXMixer: {e}")
            results['lxmixer'] = float('inf')

        print(f"Franz: {results['franz']}, LXMixer: {results['lxmixer']}")
        return results['franz'], results['lxmixer']


def saveResult(res):
    global cnots_franz, cnots_lxmixer
    cnots_franz.append(res[0])
    cnots_lxmixer.append(res[1])


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

    # Test the worker directly first
    print("Testing worker.get_costs directly...")
    try:
        test_B_strings, test_B_integers = worker.sample_B(2)
        print(f"Test B strings: {test_B_strings}")
        print(f"Test B integers: {test_B_integers}")
        test_result = worker.get_costs(test_B_strings, test_B_integers)
        print(f"Test result: {test_result}")
    except Exception as e:
        print(f"Error in direct test: {e}")
        import traceback
        traceback.print_exc()
        return

    # m = |B|
    m_list = list(range(2, 2**n + 1))

    # Arrays for Franz mixer
    min_cnots_franz = np.zeros(len(m_list))
    max_cnots_franz = np.zeros(len(m_list))
    mean_cnots_franz = np.zeros(len(m_list))
    var_cnots_franz = np.zeros(len(m_list))

    # Arrays for LXMixer
    min_cnots_lxmixer = np.zeros(len(m_list))
    max_cnots_lxmixer = np.zeros(len(m_list))
    mean_cnots_lxmixer = np.zeros(len(m_list))
    var_cnots_lxmixer = np.zeros(len(m_list))

    for i, m in enumerate(m_list):
        print(f"\nProcessing m={m} (|B|={m}) - {i+1}/{len(m_list)}")
        global cnots_franz, cnots_lxmixer
        cnots_franz = []
        cnots_lxmixer = []
        
        pool = Pool()
        results = []
        for j in range(num_samples):
            B_strings, B_integers = worker.sample_B(m)
            result = pool.apply_async(
                worker.get_costs, args=(B_strings, B_integers), callback=saveResult
            )
            results.append(result)
        pool.close()
        pool.join()
        
        # Check for exceptions in the results
        failed_count = 0
        for j, result in enumerate(results):
            try:
                result.get()  # This will raise any exception that occurred
            except Exception as e:
                print(f"Exception in worker {j}: {e}")
                failed_count += 1
                import traceback
                traceback.print_exc()
        
        print(f"Collected {len(cnots_franz)} Franz results and {len(cnots_lxmixer)} LXMixer results for m={m}")
        print(f"Failed: {failed_count}/{num_samples}")
        
        if len(cnots_franz) == 0 or len(cnots_lxmixer) == 0:
            print(f"No results collected for m={m} - all tasks failed!")
            # Set to NaN instead of skipping to maintain array indices
            min_cnots_franz[i] = np.nan
            max_cnots_franz[i] = np.nan
            mean_cnots_franz[i] = np.nan
            var_cnots_franz[i] = np.nan
            min_cnots_lxmixer[i] = np.nan
            max_cnots_lxmixer[i] = np.nan
            mean_cnots_lxmixer[i] = np.nan
            var_cnots_lxmixer[i] = np.nan
            continue

        # Calculate statistics for Franz mixer
        min_cnots_franz[i] = np.min(cnots_franz)
        max_cnots_franz[i] = np.max(cnots_franz)
        mean_cnots_franz[i] = np.mean(cnots_franz)
        var_cnots_franz[i] = np.var(cnots_franz)

        # Calculate statistics for LXMixer
        min_cnots_lxmixer[i] = np.min(cnots_lxmixer)
        max_cnots_lxmixer[i] = np.max(cnots_lxmixer)
        mean_cnots_lxmixer[i] = np.mean(cnots_lxmixer)
        var_cnots_lxmixer[i] = np.var(cnots_lxmixer)

        print(int(100 * (i + 1) / len(m_list)), "%")

    # Print final data summary
    print(f"\nFinal data summary:")
    print(f"m_list: {m_list}")
    print(f"Franz means: {mean_cnots_franz}")
    print(f"LXMixer means: {mean_cnots_lxmixer}")
    print(f"Franz has {np.sum(~np.isnan(mean_cnots_franz))} valid values")
    print(f"LXMixer has {np.sum(~np.isnan(mean_cnots_lxmixer))} valid values")

    # Create comparison plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Filter out NaN values for plotting
    valid_indices = ~(np.isnan(mean_cnots_franz) | np.isnan(mean_cnots_lxmixer))
    m_list_valid = np.array(m_list)[valid_indices]
    
    if len(m_list_valid) == 0:
        print("No valid data points to plot!")
        return
    
    mean_franz_valid = mean_cnots_franz[valid_indices]
    var_franz_valid = var_cnots_franz[valid_indices]
    min_franz_valid = min_cnots_franz[valid_indices]
    max_franz_valid = max_cnots_franz[valid_indices]
    
    mean_lxmixer_valid = mean_cnots_lxmixer[valid_indices]
    var_lxmixer_valid = var_cnots_lxmixer[valid_indices]
    min_lxmixer_valid = min_cnots_lxmixer[valid_indices]
    max_lxmixer_valid = max_cnots_lxmixer[valid_indices]

    # Plot Franz mixer results
    plt.plot(m_list_valid, mean_franz_valid, "o-b", label="Franz Mixer (mean)")
    plt.fill_between(
        m_list_valid,
        mean_franz_valid - np.sqrt(var_franz_valid),
        mean_franz_valid + np.sqrt(var_franz_valid),
        color="blue",
        alpha=0.35,
        label="Franz Mixer (1σ)"
    )
    plt.plot(m_list_valid, min_franz_valid, ":k", alpha=0.7)
    plt.plot(m_list_valid, max_franz_valid, ":k", alpha=0.7)

    # Plot LXMixer results
    plt.plot(m_list_valid, mean_lxmixer_valid, "x-r", label="LXMixer (mean)")
    plt.fill_between(
        m_list_valid,
        mean_lxmixer_valid - np.sqrt(var_lxmixer_valid),
        mean_lxmixer_valid + np.sqrt(var_lxmixer_valid),
        color="red",
        alpha=0.35,
        label="LXMixer (1σ)"
    )
    plt.plot(m_list_valid, min_lxmixer_valid, ":k", alpha=0.7)
    plt.plot(m_list_valid, max_lxmixer_valid, ":k", alpha=0.7)

    plt.legend()
    plt.xlim(min(m_list_valid), max(m_list_valid))
    plt.xlabel(r"$|B|$")
    plt.ylabel("#CNOTS")
    plt.title(f"Comparison of Franz Mixer vs LXMixer (n={n})")
    plt.grid()
    plt.savefig(f"comparison_n{n}.png")
    plt.clf()

    
if __name__ == "__main__":
    pass
    # n=dimension of Hilbert space
    for n in [3]:
        main(n)