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
        timing_results = {}
        
        # Franz mixer (uses strings) - following the exact pattern from random_B_franz_string.py
        try:
            start_time = time.time()
            
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
            
            end_time = time.time()
            
            # Use the "optimal" case (not chain, not reduced)
            results['franz'] = cnots
            timing_results['franz'] = end_time - start_time
            print(f"Franz mixer - optimal cost: {cnots}, time: {timing_results['franz']:.4f}s")
            
        except Exception as e:
            print(f"Error with Franz mixer: {e}")
            import traceback
            traceback.print_exc()
            results['franz'] = float('inf')
            timing_results['franz'] = float('inf')
        
        # LXMixer (uses integers)
        try:
            start_time = time.time()
            # if len(B_integers) == 7:
            #     print("These are the B: ",B_integers)
            mixer = LXMixer(B_integers, nL)
            mixer.compute_family_of_valid_graphs()
            mixer.compute_all_orbits()
            stabilizer = Stabilizer(B=mixer.B, n=mixer.nL, orbit_dictionary=mixer.orbits)
            stabilizer.compute_minimal_generating_sets()
            stabilizer.compute_projector_stabilizers()
            mixer.compute_costs()
            
            #best_Xs, best_Zs, best_cost = mixer.find_best_mixer()
            mixer.find_best_mixer()
            best_cost = mixer.best_cost
            
            end_time = time.time()
            
            results['lxmixer'] = best_cost
            timing_results['lxmixer'] = end_time - start_time
            print(f"LXMixer - optimal cost: {best_cost}, time: {timing_results['lxmixer']:.4f}s")
            
        except Exception as e:
            print(f"Error with LXMixer: {e}")
            # Store the failed B set for later analysis
            failed_B_info = {
                'B_strings': B_strings,
                'B_integers': B_integers,
                'error': str(e)
            }
            results['lxmixer'] = float('inf')
            timing_results['lxmixer'] = float('inf')
            # Return the failed B info as additional data
            return results['franz'], results['lxmixer'], timing_results['franz'], timing_results['lxmixer'], failed_B_info

        print(f"Franz: {results['franz']}, LXMixer: {results['lxmixer']}")
        return results['franz'], results['lxmixer'], timing_results['franz'], timing_results['lxmixer'], None


def saveResult(res):
    global cnots_franz, cnots_lxmixer, times_franz, times_lxmixer, failed_B_sets, infinity_B_sets
    cnots_franz.append(res[0])
    cnots_lxmixer.append(res[1])
    times_franz.append(res[2])
    times_lxmixer.append(res[3])
    
    # Check if there's a failed B set (5th element indicates LXMixer failure)
    if len(res) > 4 and res[4] is not None:
        failed_B_sets.append(res[4])
    
    # Check if LXMixer returned infinity (even if it didn't "fail")
    if res[1] == float('inf'):
        # We need to get the B set info - this will be handled in the main loop
        # since we don't have access to the B set here
        pass


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

    # Arrays for timing data
    min_times_franz = np.zeros(len(m_list))
    max_times_franz = np.zeros(len(m_list))
    mean_times_franz = np.zeros(len(m_list))
    var_times_franz = np.zeros(len(m_list))

    min_times_lxmixer = np.zeros(len(m_list))
    max_times_lxmixer = np.zeros(len(m_list))
    mean_times_lxmixer = np.zeros(len(m_list))
    var_times_lxmixer = np.zeros(len(m_list))

    # Track failed B sets for LXMixer
    all_failed_B_sets = []
    all_infinity_B_sets = []

    for i, m in enumerate(m_list):
        print(f"\nProcessing m={m} (|B|={m}) - {i+1}/{len(m_list)}")
        global cnots_franz, cnots_lxmixer, times_franz, times_lxmixer, failed_B_sets, infinity_B_sets
        cnots_franz = []
        cnots_lxmixer = []
        times_franz = []
        times_lxmixer = []
        failed_B_sets = []
        infinity_B_sets = []
        
        pool = Pool()
        results = []
        worker_B_sets = []  # Store B sets for each worker
        
        for j in range(num_samples):
            B_strings, B_integers = worker.sample_B(m)
            worker_B_sets.append((B_strings, B_integers))  # Store the B set
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
                
                # For failed workers, we have the actual B set that caused the failure
                failed_B_strings, failed_B_integers = worker_B_sets[j]
                print(f"Worker {j} failed with B set: {failed_B_strings}")
                
                # Create a failed B info entry
                failed_B_info = {
                    'B_strings': failed_B_strings,
                    'B_integers': failed_B_integers,
                    'error': str(e)
                }
                all_failed_B_sets.append((m, failed_B_info))
        
        print(f"Collected {len(cnots_franz)} Franz results and {len(cnots_lxmixer)} LXMixer results for m={m}")
        print(f"Failed: {failed_count}/{num_samples}")
        
        # Debug: Show what values we actually collected
        print(f"Franz results: {cnots_franz}")
        print(f"LXMixer results: {cnots_lxmixer}")
        
        # Store failed B sets for this m value
        if failed_B_sets:
            print(f"LXMixer failed on {len(failed_B_sets)} B sets for m={m}")
            all_failed_B_sets.extend([(m, failed_B) for failed_B in failed_B_sets])
        
        # Check for infinity results and store those B sets
        for j, lxmixer_result in enumerate(cnots_lxmixer):
            if lxmixer_result == float('inf'):
                # Find the corresponding B set
                if j < len(worker_B_sets):
                    B_strings, B_integers = worker_B_sets[j]
                    infinity_B_info = {
                        'B_strings': B_strings,
                        'B_integers': B_integers,
                        'error': 'LXMixer returned infinity (no exception thrown)'
                    }
                    infinity_B_sets.append(infinity_B_info)
        
        if infinity_B_sets:
            print(f"LXMixer returned infinity for {len(infinity_B_sets)} B sets for m={m}")
            all_infinity_B_sets.extend([(m, infinity_B) for infinity_B in infinity_B_sets])
        
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
            
            # Timing data
            min_times_franz[i] = np.nan
            max_times_franz[i] = np.nan
            mean_times_franz[i] = np.nan
            var_times_franz[i] = np.nan
            min_times_lxmixer[i] = np.nan
            max_times_lxmixer[i] = np.nan
            mean_times_lxmixer[i] = np.nan
            var_times_lxmixer[i] = np.nan
            continue

        # Calculate statistics for Franz mixer
        min_cnots_franz[i] = np.min(cnots_franz)
        max_cnots_franz[i] = np.max(cnots_franz)
        mean_cnots_franz[i] = np.mean(cnots_franz)
        var_cnots_franz[i] = np.var(cnots_franz)

        # Calculate statistics for LXMixer
        print(f"LXMixer results for m={m}: {cnots_lxmixer}")
        print(f"LXMixer finite values: {[x for x in cnots_lxmixer if np.isfinite(x)]}")
        
        if len([x for x in cnots_lxmixer if np.isfinite(x)]) == 0:
            print(f"All LXMixer values are infinite/NaN for m={m}")
            min_cnots_lxmixer[i] = np.nan
            max_cnots_lxmixer[i] = np.nan
            mean_cnots_lxmixer[i] = np.nan
            var_cnots_lxmixer[i] = np.nan
        else:
            min_cnots_lxmixer[i] = np.min(cnots_lxmixer)
            max_cnots_lxmixer[i] = np.max(cnots_lxmixer)
            mean_cnots_lxmixer[i] = np.mean(cnots_lxmixer)
            var_cnots_lxmixer[i] = np.var(cnots_lxmixer)

        # Calculate timing statistics for Franz mixer
        min_times_franz[i] = np.min(times_franz)
        max_times_franz[i] = np.max(times_franz)
        mean_times_franz[i] = np.mean(times_franz)
        var_times_franz[i] = np.var(times_franz)

        # Calculate timing statistics for LXMixer
        min_times_lxmixer[i] = np.min(times_lxmixer)
        max_times_lxmixer[i] = np.max(times_lxmixer)
        mean_times_lxmixer[i] = np.mean(times_lxmixer)
        var_times_lxmixer[i] = np.var(times_lxmixer)

        print(int(100 * (i + 1) / len(m_list)), "%")

    # Print final data summary
    print(f"\nFinal data summary:")
    print(f"m_list: {m_list}")
    print(f"Franz means: {mean_cnots_franz}")
    print(f"LXMixer means: {mean_cnots_lxmixer}")
    print(f"Franz time means: {mean_times_franz}")
    print(f"LXMixer time means: {mean_times_lxmixer}")
    print(f"Franz has {np.sum(~np.isnan(mean_cnots_franz))} valid values")
    print(f"LXMixer has {np.sum(~np.isnan(mean_cnots_lxmixer))} valid values")

    # Debug: Check for infinite values
    print(f"\nDebugging infinite values:")
    print(f"Franz infinite values: {np.sum(np.isinf(mean_cnots_franz))}")
    print(f"LXMixer infinite values: {np.sum(np.isinf(mean_cnots_lxmixer))}")
    
    # Check each m value individually
    for i, m in enumerate(m_list):
        franz_val = mean_cnots_franz[i]
        lxmixer_val = mean_cnots_lxmixer[i]
        print(f"m={m}: Franz={franz_val}, LXMixer={lxmixer_val}")
        if np.isnan(franz_val) or np.isnan(lxmixer_val):
            print(f"  -> NaN detected")
        if np.isinf(franz_val) or np.isinf(lxmixer_val):
            print(f"  -> Infinity detected")

    # Create comparison plot for CNOT costs and timing
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Filter out NaN and infinite values for plotting
    valid_indices = ~(np.isnan(mean_cnots_franz) | np.isnan(mean_cnots_lxmixer) | 
                      np.isinf(mean_cnots_franz) | np.isinf(mean_cnots_lxmixer))
    m_list_valid = np.array(m_list)[valid_indices]
    
    print(f"\nFiltering results:")
    print(f"Total data points: {len(m_list)}")
    print(f"Valid data points after filtering: {len(m_list_valid)}")
    print(f"Filtered out data points: {len(m_list) - len(m_list_valid)}")
    
    if len(m_list_valid) == 0:
        print("No valid data points to plot!")
        return
    
    # CNOT cost data
    mean_franz_valid = mean_cnots_franz[valid_indices]
    var_franz_valid = var_cnots_franz[valid_indices]
    min_franz_valid = min_cnots_franz[valid_indices]
    max_franz_valid = max_cnots_franz[valid_indices]
    
    mean_lxmixer_valid = mean_cnots_lxmixer[valid_indices]
    var_lxmixer_valid = var_cnots_lxmixer[valid_indices]
    min_lxmixer_valid = min_cnots_lxmixer[valid_indices]
    max_lxmixer_valid = max_cnots_lxmixer[valid_indices]

    # Timing data
    mean_times_franz_valid = mean_times_franz[valid_indices]
    var_times_franz_valid = var_times_franz[valid_indices]
    min_times_franz_valid = min_times_franz[valid_indices]
    max_times_franz_valid = max_times_franz[valid_indices]
    
    mean_times_lxmixer_valid = mean_times_lxmixer[valid_indices]
    var_times_lxmixer_valid = var_times_lxmixer[valid_indices]
    min_times_lxmixer_valid = min_times_lxmixer[valid_indices]
    max_times_lxmixer_valid = max_times_lxmixer[valid_indices]

    # Plot CNOT costs (left subplot)
    ax1.plot(m_list_valid, mean_franz_valid, "o-b", label="Franz Mixer (mean)")
    ax1.fill_between(
        m_list_valid,
        mean_franz_valid - np.sqrt(var_franz_valid),
        mean_franz_valid + np.sqrt(var_franz_valid),
        color="blue",
        alpha=0.35,
        label="Franz Mixer (1σ)"
    )
    ax1.plot(m_list_valid, min_franz_valid, ":b", alpha=0.7)
    ax1.plot(m_list_valid, max_franz_valid, ":b", alpha=0.7)

    ax1.plot(m_list_valid, mean_lxmixer_valid, "x-r", label="LXMixer (mean)")
    ax1.fill_between(
        m_list_valid,
        mean_lxmixer_valid - np.sqrt(var_lxmixer_valid),
        mean_lxmixer_valid + np.sqrt(var_lxmixer_valid),
        color="red",
        alpha=0.35,
        label="LXMixer (1σ)"
    )
    ax1.plot(m_list_valid, min_lxmixer_valid, ":r", alpha=0.7)
    ax1.plot(m_list_valid, max_lxmixer_valid, ":r", alpha=0.7)

    ax1.legend()
    ax1.set_xlim(min(m_list_valid), max(m_list_valid))
    ax1.set_xlabel(r"$|B|$")
    ax1.set_ylabel("#CNOTS")
    ax1.set_title(f"CNOT Cost Comparison (n={n})")
    ax1.grid()

    # Plot execution times (right subplot)
    ax2.plot(m_list_valid, mean_times_franz_valid, "o-b", label="Franz Mixer (mean)")
    ax2.fill_between(
        m_list_valid,
        mean_times_franz_valid - np.sqrt(var_times_franz_valid),
        mean_times_franz_valid + np.sqrt(var_times_franz_valid),
        color="blue",
        alpha=0.35,
        label="Franz Mixer (1σ)"
    )
    ax2.plot(m_list_valid, min_times_franz_valid, ":b", alpha=0.7)
    ax2.plot(m_list_valid, max_times_franz_valid, ":b", alpha=0.7)

    ax2.plot(m_list_valid, mean_times_lxmixer_valid, "x-r", label="LXMixer (mean)")
    ax2.fill_between(
        m_list_valid,
        mean_times_lxmixer_valid - np.sqrt(var_times_lxmixer_valid),
        mean_times_lxmixer_valid + np.sqrt(var_times_lxmixer_valid),
        color="red",
        alpha=0.35,
        label="LXMixer (1σ)"
    )
    ax2.plot(m_list_valid, min_times_lxmixer_valid, ":r", alpha=0.7)
    ax2.plot(m_list_valid, max_times_lxmixer_valid, ":r", alpha=0.7)

    ax2.legend()
    ax2.set_xlim(min(m_list_valid), max(m_list_valid))
    ax2.set_xlabel(r"$|B|$")
    ax2.set_ylabel("Execution Time (s)")
    ax2.set_title(f"Execution Time Comparison (n={n})")
    ax2.grid()

    plt.tight_layout()
    plt.savefig(f"comparison_n{n}.png")
    plt.clf()

    # Create separate timing plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(m_list_valid, mean_times_franz_valid, "o-b", label="Franz Mixer (mean)")
    ax.fill_between(
        m_list_valid,
        mean_times_franz_valid - np.sqrt(var_times_franz_valid),
        mean_times_franz_valid + np.sqrt(var_times_franz_valid),
        color="blue",
        alpha=0.35,
        label="Franz Mixer (1σ)"
    )
    ax.plot(m_list_valid, min_times_franz_valid, ":b", alpha=0.7)
    ax.plot(m_list_valid, max_times_franz_valid, ":b", alpha=0.7)

    ax.plot(m_list_valid, mean_times_lxmixer_valid, "x-r", label="LXMixer (mean)")
    ax.fill_between(
        m_list_valid,
        mean_times_lxmixer_valid - np.sqrt(var_times_lxmixer_valid),
        mean_times_lxmixer_valid + np.sqrt(var_times_lxmixer_valid),
        color="red",
        alpha=0.35,
        label="LXMixer (1σ)"
    )
    ax.plot(m_list_valid, min_times_lxmixer_valid, ":r", alpha=0.7)
    ax.plot(m_list_valid, max_times_lxmixer_valid, ":r", alpha=0.7)

    ax.legend()
    ax.set_xlim(min(m_list_valid), max(m_list_valid))
    ax.set_xlabel(r"$|B|$")
    ax.set_ylabel("Execution Time (s)")
    ax.set_title(f"Execution Time Comparison: Franz Mixer vs LXMixer (n={n})")
    ax.grid()
    
    plt.savefig(f"timing_comparison_n{n}.png")
    plt.clf()

    # Print failed B sets at the very end
    print("\n" + "="*80)
    print("FAILED B SETS FOR LXMixer:")
    print("="*80)
    
    if all_failed_B_sets:
        for m, failed_B_info in all_failed_B_sets:
            print(f"\nFailed B set for m={m} (|B|={m}):")
            print(f"  B_strings: {failed_B_info['B_strings']}")
            print(f"  B_integers: {failed_B_info['B_integers']}")
            print(f"  Error: {failed_B_info['error']}")
            print("-" * 40)
        
        print(f"\nTotal failed B sets: {len(all_failed_B_sets)}")
        
        # Group by m value
        m_failures = {}
        for m, _ in all_failed_B_sets:
            m_failures[m] = m_failures.get(m, 0) + 1
        
        print("Failures by |B| value:")
        for m in sorted(m_failures.keys()):
            print(f"  |B|={m}: {m_failures[m]} failures")
    else:
        print("No failed B sets found!")
    
    print("="*80)

    # Print infinity B sets
    print("\n" + "="*80)
    print("B SETS WHERE LXMixer RETURNED INFINITY:")
    print("="*80)
    
    if all_infinity_B_sets:
        for m, infinity_B_info in all_infinity_B_sets:
            print(f"\nInfinity B set for m={m} (|B|={m}):")
            print(f"  B_strings: {infinity_B_info['B_strings']}")
            print(f"  B_integers: {infinity_B_info['B_integers']}")
            print(f"  Error: {infinity_B_info['error']}")
            print("-" * 40)
        
        print(f"\nTotal infinity B sets: {len(all_infinity_B_sets)}")
        
        # Group by m value
        m_infinities = {}
        for m, _ in all_infinity_B_sets:
            m_infinities[m] = m_infinities.get(m, 0) + 1
        
        print("Infinity results by |B| value:")
        for m in sorted(m_infinities.keys()):
            print(f"  |B|={m}: {m_infinities[m]} infinity results")
    else:
        print("No infinity B sets found!")
    
    print("="*80)

    # Print combined summary
    print("\n" + "="*80)
    print("SUMMARY - ALL PROBLEMATIC B SETS:")
    print("="*80)
    total_problematic = len(all_failed_B_sets) + len(all_infinity_B_sets)
    print(f"Total problematic B sets: {total_problematic}")
    print(f"  - Failed with exceptions: {len(all_failed_B_sets)}")
    print(f"  - Returned infinity: {len(all_infinity_B_sets)}")
    print("="*80)

    
if __name__ == "__main__":
    pass
    # n=dimension of Hilbert space
    for n in [3]:
        main(n)