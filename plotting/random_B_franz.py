from multiprocessing import Pool

#import tikzplotlib

import itertools
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving plots
from matplotlib import pyplot as plt
import numpy as np
from random import sample
import time
import sys
import os
import psutil  # For system resource monitoring
import gc  # For garbage collection monitoring

# Add parent directory to path to import Mixer
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mixer_path = os.path.join(parent_dir, 'Mixer.py')
sys.path.append(parent_dir)

from Mixer import *

# Importing necessary classes from Mixer module
sys.path.append(r"C:\Users\sanne\LogicalXMixer")
from Mixer_Franz import MixerFranz

# CONFIGURATION: Choose which mixers to run
# Set to True to run that mixer, False to skip it
RUN_FRANZ_MIXER = False    # Set to False to skip Franz mixer
RUN_LXMIXER = True        # Set to False to skip LXMixer

# CONFIGURATION: Simple optimizations
USE_SIMPLE_MONITORING = True  # Keep basic timing analysis without heavy overhead

# You can also run only one mixer by setting one to False:
# RUN_FRANZ_MIXER = True; RUN_LXMIXER = False   # Only Franz
# RUN_FRANZ_MIXER = False; RUN_LXMIXER = True   # Only LXMixer

# USAGE EXAMPLES:
# 1. Run both mixers (default):
#    RUN_FRANZ_MIXER = True; RUN_LXMIXER = True
# 2. Only run Franz mixer (faster, no LXMixer computation):
#    RUN_FRANZ_MIXER = True; RUN_LXMIXER = False
# 3. Only run LXMixer (faster, no Franz mixer computation):
#    RUN_FRANZ_MIXER = False; RUN_LXMIXER = True


class Worker:

    """n = dim(hilbert space)"""

    def __init__(self, n, checkmixer=False, run_franz=True, run_lxmixer=True):
        self.n = n  # Store n for later use
        self.all_states = ["".join(i) for i in itertools.product("01", repeat=n)]
        self.checkmixer = checkmixer
        self.run_franz = run_franz
        self.run_lxmixer = run_lxmixer

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
        
        # Initialize with None values
        results['franz'] = None
        results['lxmixer'] = None
        timing_results['franz'] = 0.0
        timing_results['lxmixer'] = 0.0
        
        # Franz mixer (uses strings) - following the exact pattern from random_B_franz_string.py
        if self.run_franz:
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
        else:
            # Skip Franz mixer - use default value
            results['franz'] = 0  # or float('inf') if you prefer
            timing_results['franz'] = 0.0
        
        # LXMixer (uses integers)
        if self.run_lxmixer:
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
        else:
            # Skip LXMixer - use default value
            results['lxmixer'] = 0  # or float('inf') if you prefer
            timing_results['lxmixer'] = 0.0

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


def get_cpu_temperature():
    """Try to get CPU temperature if available"""
    try:
        temps = psutil.sensors_temperatures()
        if temps:
            # Try common temperature sensor names
            for name in ['coretemp', 'cpu_thermal', 'acpi']:
                if name in temps:
                    return temps[name][0].current
        return None
    except:
        return None

def main(n, num_samples=100):
    print("n=", n)
    print(f"Configuration: Franz Mixer={'ON' if RUN_FRANZ_MIXER else 'OFF'}, LXMixer={'ON' if RUN_LXMIXER else 'OFF'}")
    
    # Validate configuration
    if not RUN_FRANZ_MIXER and not RUN_LXMIXER:
        print("ERROR: Both mixers are disabled! Please enable at least one mixer.")
        return
    
    worker = Worker(n, checkmixer=False, run_franz=RUN_FRANZ_MIXER, run_lxmixer=RUN_LXMIXER)

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
    
    # # RESUME: Skip completed m values - change this number to resume from specific m
    # start_from_m = 11  # CHANGE THIS: Set to the m value you want to resume from
    # try:
    #     start_index = m_list.index(start_from_m)
    #     print(f"RESUMING from m={start_from_m} (skipping {start_index} completed values)")
    #     m_list = m_list[start_index:]
    # except ValueError:
    #     print(f"Error: m={start_from_m} not found in m_list")
    #     return

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
        
        # Lightweight resource monitoring (only if enabled)
        if USE_SIMPLE_MONITORING:
            process = psutil.Process()
            memory_before = process.memory_info().rss / 1024 / 1024  # MB
            print(f"Process memory: {memory_before:.1f} MB")
        
        # Force garbage collection before starting
        gc.collect()
        
        global cnots_franz, cnots_lxmixer, times_franz, times_lxmixer, failed_B_sets, infinity_B_sets
        cnots_franz = []
        cnots_lxmixer = []
        times_franz = []
        times_lxmixer = []
        failed_B_sets = []
        infinity_B_sets = []
        
        # Time the entire batch processing
        batch_start_time = time.time()
        
        # Simple single pool processing (back to original, efficient method)
        print(f"Processing all {num_samples} workers in single pool")
        pool = Pool()
        results = []
        worker_B_sets = []  # Store B sets for each worker
        
        # Track individual worker timings (only if monitoring enabled)
        worker_start_times = []
        
        for j in range(num_samples):
            if USE_SIMPLE_MONITORING:
                worker_start = time.time()
                worker_start_times.append(worker_start)
            
            B_strings, B_integers = worker.sample_B(m)
            worker_B_sets.append((B_strings, B_integers))  # Store the B set
            result = pool.apply_async(
                worker.get_costs, args=(B_strings, B_integers), callback=saveResult
            )
            results.append(result)
        
        pool.close()
        
        print(f"All {num_samples} workers started. Waiting for completion...")
        join_start = time.time()
        pool.join()
        join_end = time.time()
        
        batch_end_time = time.time()
        total_batch_time = batch_end_time - batch_start_time
        join_time = join_end - join_start
        
        print(f"Batch completed in {total_batch_time:.1f}s (join took {join_time:.1f}s)")
        
        # Check system resources after completion (lightweight)
        if USE_SIMPLE_MONITORING:
            memory_after = process.memory_info().rss / 1024 / 1024  # MB
            print(f"Process memory after: {memory_after:.1f} MB (change: {memory_after - memory_before:+.1f} MB)")
        
        # Check for exceptions in the results and track completion times
        failed_count = 0
        worker_completion_times = []
        
        for j, result in enumerate(results):
            try:
                result.get()  # This will raise any exception that occurred
                
                # Only track timing if monitoring is enabled and we have start times
                if USE_SIMPLE_MONITORING and j < len(worker_start_times):
                    completion_time = time.time() - worker_start_times[j]
                    worker_completion_times.append(completion_time)
                    
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
        
        # Analyze worker timing patterns (only if monitoring enabled and we have data)
        if USE_SIMPLE_MONITORING and worker_completion_times and len(worker_completion_times) >= 10:
            avg_time = np.mean(worker_completion_times)
            std_time = np.std(worker_completion_times)
            
            # Simple analysis: first vs last workers
            if len(worker_completion_times) >= 20:
                first_10_avg = np.mean(worker_completion_times[:10])
                last_10_avg = np.mean(worker_completion_times[-10:])
                
                print(f"Worker timing analysis:")
                print(f"  First 10 workers avg: {first_10_avg:.1f}s")
                print(f"  Last 10 workers avg: {last_10_avg:.1f}s")
                print(f"  Overall avg: {avg_time:.1f}s ± {std_time:.1f}s")
                
                # Simple interpretation
                if last_10_avg > first_10_avg * 1.5:
                    print(f"  ⚠️  WARNING: Resource exhaustion detected (last workers {last_10_avg/first_10_avg:.1f}x slower)")
                elif first_10_avg > last_10_avg * 1.3:
                    print(f"  ✅ INFO: Normal startup effect (first workers {first_10_avg/last_10_avg:.1f}x slower)")
                else:
                    print(f"  ✅ GOOD: Consistent worker timing")
        
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
    
    if RUN_FRANZ_MIXER:
        print(f"Franz means: {mean_cnots_franz}")
        print(f"Franz time means: {mean_times_franz}")
        print(f"Franz has {np.sum(~np.isnan(mean_cnots_franz))} valid values")
    else:
        print("Franz mixer was disabled")
    
    if RUN_LXMIXER:
        print(f"LXMixer means: {mean_cnots_lxmixer}")
        print(f"LXMixer time means: {mean_times_lxmixer}")
        print(f"LXMixer has {np.sum(~np.isnan(mean_cnots_lxmixer))} valid values")
    else:
        print("LXMixer was disabled")

    # Debug: Check for infinite values
    print(f"\nDebugging infinite values:")
    if RUN_FRANZ_MIXER:
        print(f"Franz infinite values: {np.sum(np.isinf(mean_cnots_franz))}")
    if RUN_LXMIXER:
        print(f"LXMixer infinite values: {np.sum(np.isinf(mean_cnots_lxmixer))}")
    
    # Check each m value individually
    for i, m in enumerate(m_list):
        if RUN_FRANZ_MIXER and RUN_LXMIXER:
            franz_val = mean_cnots_franz[i]
            lxmixer_val = mean_cnots_lxmixer[i]
            print(f"m={m}: Franz={franz_val}, LXMixer={lxmixer_val}")
            if np.isnan(franz_val) or np.isnan(lxmixer_val):
                print(f"  -> NaN detected")
            if np.isinf(franz_val) or np.isinf(lxmixer_val):
                print(f"  -> Infinity detected")
        elif RUN_FRANZ_MIXER:
            franz_val = mean_cnots_franz[i]
            print(f"m={m}: Franz={franz_val}")
            if np.isnan(franz_val):
                print(f"  -> NaN detected")
            if np.isinf(franz_val):
                print(f"  -> Infinity detected")
        elif RUN_LXMIXER:
            lxmixer_val = mean_cnots_lxmixer[i]
            print(f"m={m}: LXMixer={lxmixer_val}")
            if np.isnan(lxmixer_val):
                print(f"  -> NaN detected")
            if np.isinf(lxmixer_val):
                print(f"  -> Infinity detected")

    # Create comparison plot for CNOT costs and timing
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Filter out NaN and infinite values for plotting based on enabled mixers
    if RUN_FRANZ_MIXER and RUN_LXMIXER:
        # Both mixers enabled - filter based on both
        valid_indices = ~(np.isnan(mean_cnots_franz) | np.isnan(mean_cnots_lxmixer) | 
                          np.isinf(mean_cnots_franz) | np.isinf(mean_cnots_lxmixer))
    elif RUN_FRANZ_MIXER:
        # Only Franz mixer enabled
        valid_indices = ~(np.isnan(mean_cnots_franz) | np.isinf(mean_cnots_franz))
    elif RUN_LXMIXER:
        # Only LXMixer enabled
        valid_indices = ~(np.isnan(mean_cnots_lxmixer) | np.isinf(mean_cnots_lxmixer))
    else:
        print("No mixers enabled - cannot create plots!")
        return
    
    m_list_valid = np.array(m_list)[valid_indices]
    
    print(f"\nFiltering results:")
    print(f"Total data points: {len(m_list)}")
    print(f"Valid data points after filtering: {len(m_list_valid)}")
    print(f"Filtered out data points: {len(m_list) - len(m_list_valid)}")
    
    if len(m_list_valid) == 0:
        print("No valid data points to plot!")
        return
    
    # CNOT cost data
    if RUN_FRANZ_MIXER:
        mean_franz_valid = mean_cnots_franz[valid_indices]
        var_franz_valid = var_cnots_franz[valid_indices]
        min_franz_valid = min_cnots_franz[valid_indices]
        max_franz_valid = max_cnots_franz[valid_indices]
    
    if RUN_LXMIXER:
        mean_lxmixer_valid = mean_cnots_lxmixer[valid_indices]
        var_lxmixer_valid = var_cnots_lxmixer[valid_indices]
        min_lxmixer_valid = min_cnots_lxmixer[valid_indices]
        max_lxmixer_valid = max_cnots_lxmixer[valid_indices]

    # Timing data
    if RUN_FRANZ_MIXER:
        mean_times_franz_valid = mean_times_franz[valid_indices]
        var_times_franz_valid = var_times_franz[valid_indices]
        min_times_franz_valid = min_times_franz[valid_indices]
        max_times_franz_valid = max_times_franz[valid_indices]
    
    if RUN_LXMIXER:
        mean_times_lxmixer_valid = mean_times_lxmixer[valid_indices]
        var_times_lxmixer_valid = var_times_lxmixer[valid_indices]
        min_times_lxmixer_valid = min_times_lxmixer[valid_indices]
        max_times_lxmixer_valid = max_times_lxmixer[valid_indices]

    # Plot CNOT costs (left subplot)
    if RUN_FRANZ_MIXER:
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

    if RUN_LXMIXER:
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
    
    # Update title based on enabled mixers
    if RUN_FRANZ_MIXER and RUN_LXMIXER:
        ax1.set_title(f"CNOT Cost Comparison (n={n})")
    elif RUN_FRANZ_MIXER:
        ax1.set_title(f"Franz Mixer CNOT Cost (n={n})")
    else:
        ax1.set_title(f"LXMixer CNOT Cost (n={n})")
    ax1.grid()

    # Plot execution times (right subplot)
    if RUN_FRANZ_MIXER:
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

    if RUN_LXMIXER:
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
    
    # Update title based on enabled mixers
    if RUN_FRANZ_MIXER and RUN_LXMIXER:
        ax2.set_title(f"Execution Time Comparison (n={n})")
    elif RUN_FRANZ_MIXER:
        ax2.set_title(f"Franz Mixer Execution Time (n={n})")
    else:
        ax2.set_title(f"LXMixer Execution Time (n={n})")
    ax2.grid()

    plt.tight_layout()
    
    # Update filename based on enabled mixers
    if RUN_FRANZ_MIXER and RUN_LXMIXER:
        plt.savefig(f"comparison_n{n}.png")
    elif RUN_FRANZ_MIXER:
        plt.savefig(f"franz_only_n{n}.png")
    else:
        plt.savefig(f"lxmixer_only_n{n}.png")
    plt.clf()

    # Create separate timing plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if RUN_FRANZ_MIXER:
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

    if RUN_LXMIXER:
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
    
    # Update title based on enabled mixers
    if RUN_FRANZ_MIXER and RUN_LXMIXER:
        ax.set_title(f"Execution Time Comparison: Franz Mixer vs LXMixer (n={n})")
        timing_filename = f"timing_comparison_n{n}.png"
    elif RUN_FRANZ_MIXER:
        ax.set_title(f"Franz Mixer Execution Time (n={n})")
        timing_filename = f"timing_franz_only_n{n}.png"
    else:
        ax.set_title(f"LXMixer Execution Time (n={n})")
        timing_filename = f"timing_lxmixer_only_n{n}.png"
    
    ax.grid()
    
    plt.savefig(timing_filename)
    plt.clf()

    # Print failed B sets at the very end (only if LXMixer was enabled)
    if RUN_LXMIXER:
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
                print(f"  Note: {infinity_B_info['error']}")
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
    else:
        print("\n" + "="*80)
        print("LXMixer was disabled - no failure analysis to report")
        print("="*80)

    
if __name__ == "__main__":
    # n=dimension of Hilbert space
    
    # Examples of how to use the configuration:
    # 
    # 1. Run both mixers (default configuration above):
    #    for n in [3]:
    #        main(n)
    #
    # 2. Run only Franz mixer (modify RUN_LXMIXER = False above):
    #    for n in [3]:
    #        main(n)
    #
    # 3. Run only LXMixer (modify RUN_FRANZ_MIXER = False above):
    #    for n in [3]:
    #        main(n)
    
    for n in [4]:
        main(n)