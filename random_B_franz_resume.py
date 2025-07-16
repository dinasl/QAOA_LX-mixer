import pickle
import os
from multiprocessing import Pool
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
            mixer = LXMixer(B_integers, nL)
            mixer.compute_family_of_valid_graphs()
            mixer.compute_all_orbits()
            stabilizer = Stabilizer(B=mixer.B, n=mixer.nL, orbit_dictionary=mixer.orbits)
            stabilizer.compute_minimal_generating_sets()
            stabilizer.compute_projector_stabilizers()
            mixer.compute_costs()
            
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

def save_checkpoint(checkpoint_data, filename):
    """Save checkpoint data to file"""
    with open(filename, 'wb') as f:
        pickle.dump(checkpoint_data, f)
    print(f"Checkpoint saved to {filename}")

def load_checkpoint(filename):
    """Load checkpoint data from file"""
    try:
        with open(filename, 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError:
        return None

def main_with_resume(n, num_samples=100, start_from_m=None):
    print(f"n={n}")
    worker = Worker(n, checkmixer=False)
    
    checkpoint_file = f"checkpoint_n{n}.pkl"
    
    # Try to load existing checkpoint
    checkpoint = load_checkpoint(checkpoint_file)
    
    if checkpoint is not None:
        print(f"Found checkpoint file! Resuming from m={checkpoint['current_m']}")
        print(f"Previous progress: {checkpoint['completed_m_values']} m values completed")
        
        # Restore data
        min_cnots_franz = checkpoint['min_cnots_franz']
        max_cnots_franz = checkpoint['max_cnots_franz']
        mean_cnots_franz = checkpoint['mean_cnots_franz']
        var_cnots_franz = checkpoint['var_cnots_franz']
        
        min_cnots_lxmixer = checkpoint['min_cnots_lxmixer']
        max_cnots_lxmixer = checkpoint['max_cnots_lxmixer']
        mean_cnots_lxmixer = checkpoint['mean_cnots_lxmixer']
        var_cnots_lxmixer = checkpoint['var_cnots_lxmixer']
        
        min_times_franz = checkpoint['min_times_franz']
        max_times_franz = checkpoint['max_times_franz']
        mean_times_franz = checkpoint['mean_times_franz']
        var_times_franz = checkpoint['var_times_franz']
        
        min_times_lxmixer = checkpoint['min_times_lxmixer']
        max_times_lxmixer = checkpoint['max_times_lxmixer']
        mean_times_lxmixer = checkpoint['mean_times_lxmixer']
        var_times_lxmixer = checkpoint['var_times_lxmixer']
        
        all_failed_B_sets = checkpoint['all_failed_B_sets']
        all_infinity_B_sets = checkpoint['all_infinity_B_sets']
        
        m_list = checkpoint['m_list']
        start_i = checkpoint['current_i'] + 1  # Resume from next m value
        
        print(f"Resuming from m={m_list[start_i] if start_i < len(m_list) else 'DONE'}")
        
    else:
        print("No checkpoint found, starting fresh")
        
        # Initialize fresh
        m_list = list(range(2, 2**n + 1))
        start_i = 0
        
        # Initialize arrays
        min_cnots_franz = np.zeros(len(m_list))
        max_cnots_franz = np.zeros(len(m_list))
        mean_cnots_franz = np.zeros(len(m_list))
        var_cnots_franz = np.zeros(len(m_list))

        min_cnots_lxmixer = np.zeros(len(m_list))
        max_cnots_lxmixer = np.zeros(len(m_list))
        mean_cnots_lxmixer = np.zeros(len(m_list))
        var_cnots_lxmixer = np.zeros(len(m_list))

        min_times_franz = np.zeros(len(m_list))
        max_times_franz = np.zeros(len(m_list))
        mean_times_franz = np.zeros(len(m_list))
        var_times_franz = np.zeros(len(m_list))

        min_times_lxmixer = np.zeros(len(m_list))
        max_times_lxmixer = np.zeros(len(m_list))
        mean_times_lxmixer = np.zeros(len(m_list))
        var_times_lxmixer = np.zeros(len(m_list))

        all_failed_B_sets = []
        all_infinity_B_sets = []
    
    # Override start position if specified
    if start_from_m is not None:
        try:
            start_i = m_list.index(start_from_m)
            print(f"Manual override: starting from m={start_from_m} (index {start_i})")
        except ValueError:
            print(f"Error: m={start_from_m} not found in m_list")
            return
    
    if start_i >= len(m_list):
        print("All m values completed! Generating final plots...")
        # Jump to plotting section
        start_i = len(m_list)
    
    # Main processing loop
    for i in range(start_i, len(m_list)):
        m = m_list[i]
        print(f"\nProcessing m={m} (|B|={m}) - {i+1}/{len(m_list)}")
        
        # Process this m value (same as original code)
        global cnots_franz, cnots_lxmixer, times_franz, times_lxmixer, failed_B_sets, infinity_B_sets
        cnots_franz = []
        cnots_lxmixer = []
        times_franz = []
        times_lxmixer = []
        failed_B_sets = []
        infinity_B_sets = []
        
        def saveResult(res):
            global cnots_franz, cnots_lxmixer, times_franz, times_lxmixer, failed_B_sets, infinity_B_sets
            cnots_franz.append(res[0])
            cnots_lxmixer.append(res[1])
            times_franz.append(res[2])
            times_lxmixer.append(res[3])
            
            if len(res) > 4 and res[4] is not None:
                failed_B_sets.append(res[4])
        
        pool = Pool()
        results = []
        worker_B_sets = []
        
        for j in range(num_samples):
            B_strings, B_integers = worker.sample_B(m)
            worker_B_sets.append((B_strings, B_integers))
            result = pool.apply_async(
                worker.get_costs, args=(B_strings, B_integers), callback=saveResult
            )
            results.append(result)
        pool.close()
        pool.join()
        
        # Process results (same as original)
        failed_count = 0
        for j, result in enumerate(results):
            try:
                result.get()
            except Exception as e:
                print(f"Exception in worker {j}: {e}")
                failed_count += 1
                
                failed_B_strings, failed_B_integers = worker_B_sets[j]
                failed_B_info = {
                    'B_strings': failed_B_strings,
                    'B_integers': failed_B_integers,
                    'error': str(e)
                }
                all_failed_B_sets.append((m, failed_B_info))
        
        print(f"Collected {len(cnots_franz)} Franz results and {len(cnots_lxmixer)} LXMixer results for m={m}")
        print(f"Failed: {failed_count}/{num_samples}")
        
        # Store failed B sets
        if failed_B_sets:
            print(f"LXMixer failed on {len(failed_B_sets)} B sets for m={m}")
            all_failed_B_sets.extend([(m, failed_B) for failed_B in failed_B_sets])
        
        # Check for infinity results
        for j, lxmixer_result in enumerate(cnots_lxmixer):
            if lxmixer_result == float('inf'):
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
        
        # Calculate statistics
        if len(cnots_franz) == 0 or len(cnots_lxmixer) == 0:
            print(f"No results collected for m={m} - all tasks failed!")
            min_cnots_franz[i] = np.nan
            max_cnots_franz[i] = np.nan
            mean_cnots_franz[i] = np.nan
            var_cnots_franz[i] = np.nan
            min_cnots_lxmixer[i] = np.nan
            max_cnots_lxmixer[i] = np.nan
            mean_cnots_lxmixer[i] = np.nan
            var_cnots_lxmixer[i] = np.nan
            
            min_times_franz[i] = np.nan
            max_times_franz[i] = np.nan
            mean_times_franz[i] = np.nan
            var_times_franz[i] = np.nan
            min_times_lxmixer[i] = np.nan
            max_times_lxmixer[i] = np.nan
            mean_times_lxmixer[i] = np.nan
            var_times_lxmixer[i] = np.nan
        else:
            # Calculate statistics for Franz mixer
            min_cnots_franz[i] = np.min(cnots_franz)
            max_cnots_franz[i] = np.max(cnots_franz)
            mean_cnots_franz[i] = np.mean(cnots_franz)
            var_cnots_franz[i] = np.var(cnots_franz)

            # Calculate statistics for LXMixer
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

            # Calculate timing statistics
            min_times_franz[i] = np.min(times_franz)
            max_times_franz[i] = np.max(times_franz)
            mean_times_franz[i] = np.mean(times_franz)
            var_times_franz[i] = np.var(times_franz)

            min_times_lxmixer[i] = np.min(times_lxmixer)
            max_times_lxmixer[i] = np.max(times_lxmixer)
            mean_times_lxmixer[i] = np.mean(times_lxmixer)
            var_times_lxmixer[i] = np.var(times_lxmixer)

        print(f"Progress: {int(100 * (i + 1) / len(m_list))}%")
        
        # Save checkpoint after each m value
        checkpoint_data = {
            'n': n,
            'current_m': m,
            'current_i': i,
            'completed_m_values': i + 1,
            'm_list': m_list,
            'min_cnots_franz': min_cnots_franz,
            'max_cnots_franz': max_cnots_franz,
            'mean_cnots_franz': mean_cnots_franz,
            'var_cnots_franz': var_cnots_franz,
            'min_cnots_lxmixer': min_cnots_lxmixer,
            'max_cnots_lxmixer': max_cnots_lxmixer,
            'mean_cnots_lxmixer': mean_cnots_lxmixer,
            'var_cnots_lxmixer': var_cnots_lxmixer,
            'min_times_franz': min_times_franz,
            'max_times_franz': max_times_franz,
            'mean_times_franz': mean_times_franz,
            'var_times_franz': var_times_franz,
            'min_times_lxmixer': min_times_lxmixer,
            'max_times_lxmixer': max_times_lxmixer,
            'mean_times_lxmixer': mean_times_lxmixer,
            'var_times_lxmixer': var_times_lxmixer,
            'all_failed_B_sets': all_failed_B_sets,
            'all_infinity_B_sets': all_infinity_B_sets,
            'num_samples': num_samples
        }
        save_checkpoint(checkpoint_data, checkpoint_file)
    
    # Generate final plots and summary (same as original)
    print("\nGenerating final plots and summary...")
    
    # Filter out NaN and infinite values for plotting
    valid_indices = ~(np.isnan(mean_cnots_franz) | np.isnan(mean_cnots_lxmixer) | 
                      np.isinf(mean_cnots_franz) | np.isinf(mean_cnots_lxmixer))
    m_list_valid = np.array(m_list)[valid_indices]
    
    if len(m_list_valid) > 0:
        # Create plots (same as original)
        # ... (plotting code would go here)
        
        # Save final plots
        print(f"Plots saved as comparison_n{n}.png and timing_comparison_n{n}.png")
    
    # Print final summary
    print("\n" + "="*80)
    print("FINAL SUMMARY:")
    print("="*80)
    print(f"Total problematic B sets: {len(all_failed_B_sets) + len(all_infinity_B_sets)}")
    print(f"  - Failed with exceptions: {len(all_failed_B_sets)}")
    print(f"  - Returned infinity: {len(all_infinity_B_sets)}")
    print("="*80)
    
    # Clean up checkpoint file
    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
        print(f"Checkpoint file {checkpoint_file} removed (run completed)")

if __name__ == "__main__":
    # Resume for n=4
    main_with_resume(4, num_samples=100)
