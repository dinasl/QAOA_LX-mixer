from multiprocessing import Pool
import itertools
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving plots
from matplotlib import pyplot as plt
import numpy as np
from random import sample
import time

from pathlib import Path
import urllib.request
import tempfile
import os
import zipfile
import shutil

def download_github_repository():
    """Download the entire GitHub repository as a ZIP file and extract it, with caching"""
    
    # Create cache directory if it doesn't exist
    os.makedirs(CACHE_DIR, exist_ok=True)
    cached_repo_path = os.path.join(CACHE_DIR, "LogicalXMixer")
    cache_timestamp_file = os.path.join(CACHE_DIR, "last_download.txt")
    
    # Check if we have a cached version and if it's still fresh
    if os.path.exists(cached_repo_path) and os.path.exists(cache_timestamp_file):
        try:
            with open(cache_timestamp_file, 'r') as f:
                last_download_time = float(f.read().strip())
            
            current_time = time.time()
            hours_since_download = (current_time - last_download_time) / 3600
            
            if hours_since_download < CACHE_EXPIRY_HOURS:
                return cached_repo_path
            else:
                print(f"Cached repository is {hours_since_download:.1f} hours old, re-downloading...")
        except (FileNotFoundError, ValueError) as e:
            print(f"Error reading cache timestamp, re-downloading... ({e})")
    
    # print("Downloading entire Original LXMixer repository from GitHub...")
    
    # Create a temp directory for the download
    temp_dir = tempfile.mkdtemp()
    zip_path = os.path.join(temp_dir, "repo.zip")
    
    try:
        # Download the repository ZIP file
        urllib.request.urlretrieve(GITHUB_REPO_ZIP_URL, zip_path)
        
        # Extract the ZIP file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
        
        # Find the extracted repository folder (usually has a name like "LogicalXMixer-main")
        extracted_folders = [f for f in os.listdir(temp_dir) if os.path.isdir(os.path.join(temp_dir, f))]
        if not extracted_folders:
            print("Error: No folders found in extracted ZIP")
            return None
        
        repo_folder = os.path.join(temp_dir, extracted_folders[0])
        
        # Remove old cached version if it exists
        if os.path.exists(cached_repo_path):
            shutil.rmtree(cached_repo_path)
        
        # Move the downloaded repository to the cache location
        shutil.move(repo_folder, cached_repo_path)
        
        # Save the download timestamp
        with open(cache_timestamp_file, 'w') as f:
            f.write(str(time.time()))
        
        # print(f"  Repository cached to: {cached_repo_path}")
        # print(f"  Files in repository: {os.listdir(cached_repo_path)}")
        
        # Clean up the temporary directory
        shutil.rmtree(temp_dir)
        
        return cached_repo_path
        
    except Exception as e:
        print(f"Error downloading/extracting repository: {e}")
        # Clean up temp directory on error
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        return None

def force_refresh_cache():
    """Force refresh the cached repository by removing it"""
    cache_timestamp_file = os.path.join(CACHE_DIR, "last_download.txt")
    if os.path.exists(cache_timestamp_file):
        os.remove(cache_timestamp_file)
        print("Cache invalidated - next run will download fresh repository")
    else:
        print("No cache found to refresh")

def load_module(alias_name, file_path):
    import importlib.util
    import sys
    import os
    
    # Add the directory containing the module to sys.path temporarily
    module_dir = os.path.dirname(str(file_path))
    if module_dir not in sys.path:
        sys.path.insert(0, module_dir)
    
    spec = importlib.util.spec_from_file_location(alias_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias_name] = module
    spec.loader.exec_module(module)
    return module

# CONFIGURATION: Choose which mixers to run
# Set to True to run that mixer, False to skip it
RUN_ORIGINAL_LXMIXER = True    
RUN_LXMIXER_LARGEST_ORBIT = True        
RUN_LXMIXER_ALL_SUBORBIT = True         
RUN_LXMIXER_SEMI_RESTRICTED_SUBORBIT = False  # DO NOT SET TO True, not implemented yet

# Helper variable for backward compatibility
RUN_LXMIXER = RUN_LXMIXER_LARGEST_ORBIT or RUN_LXMIXER_ALL_SUBORBIT or RUN_LXMIXER_SEMI_RESTRICTED_SUBORBIT

# Get the directory where the current file is located
current_path = Path(__file__).resolve().parent

# Go one level up
mixer_new_path = current_path.parent

# GitHub configuration for mixer_old
GITHUB_REPO_ZIP_URL = "https://github.com/OpenQuantumComputing/LogicalXMixer/archive/refs/heads/main.zip"
CACHE_DIR = os.path.join(os.path.expanduser("~"), ".qaoa_mixer_cache")
CACHE_EXPIRY_HOURS = 24  # Re-download after 24 hours

# CONFIGURATION: Cache settings
FORCE_REFRESH_CACHE = False  # Set to True to force download fresh repository (ignores cache)

# Download and load mixer_old from GitHub repository
if FORCE_REFRESH_CACHE:
    force_refresh_cache()

repo_folder = download_github_repository()
if repo_folder is None:
    print("Failed to download Original LXMixer repository from GitHub!")
    exit(1)

mixer_old_path = os.path.join(repo_folder, "Mixer.py")
if not os.path.exists(mixer_old_path):
    print("Mixer.py not found in downloaded repository!")
    print(f"Available files: {os.listdir(repo_folder)}")
    exit(1)

mixer_new = load_module("Mixer", mixer_new_path / "Mixer.py")
mixer_old = load_module("Original_LXMixer", mixer_old_path)

OriginalLXMixer = mixer_old.Mixer
LXMixer = mixer_new.LXMixer



class Worker:

    """n = dim(hilbert space)"""

    def __init__(self, n, checkmixer=False, run_original_lxmixer=True, lxmixer_methods=[]):
        self.n = n  # Store n for later use
        self.all_states = ["".join(i) for i in itertools.product("01", repeat=n)]
        self.checkmixer = checkmixer
        self.run_original_lxmixer = run_original_lxmixer
        self.lxmixer_methods = lxmixer_methods  # List of LXMixer methods to run

    def sample_B(self, m):
        """m = |B|"""
        binary_strings = sample(self.all_states, m)
        binary_integers = [int(binary_str, 2) for binary_str in binary_strings]
        # Return both formats for flexibility
        return binary_strings, binary_integers
    

    def get_costs(self, B_strings, B_integers):
        nL = self.n  # Use the stored number of qubits
        
        # Finding the cnot cost for both Original LXMixer (strings) and new LXMixer (integers)
        results = {}
        timing_results = {}
        
        # Initialize with None values
        results['original_lxmixer'] = None
        results['lxmixer'] = None
        timing_results['original_lxmixer'] = 0.0
        timing_results['lxmixer'] = 0.0
        
        # Original LXMixer (uses strings)
        if self.run_original_lxmixer:
            try:
                cnots = None
                cnots_reduced = None
                
                for chain in [True, False]:
                    for reduced in [True, False]:
                        start_time = time.time()  # Time each variant separately
                        mixer_original_lxmixer = OriginalLXMixer(B_strings, reduced=reduced)
                        
                        if chain:
                            mixer_original_lxmixer.get_chain_mixer()
                        else:
                            mixer_original_lxmixer.get_best_mixer_commuting_graphs()
                        
                        end_time = time.time()
                        variant_time = end_time - start_time
                        
                        if chain and reduced:
                            cnots_chain_reduced = mixer_original_lxmixer.solution_chain_reduced_cost
                        elif chain and (not reduced):
                            cnots_chain = mixer_original_lxmixer.solution_chain_cost
                        elif (not chain) and reduced:
                            cnots_reduced = mixer_original_lxmixer.solution_reduced_cost
                        else:
                            cnots = mixer_original_lxmixer.solution_cost
                            # Store timing only for the optimal case (not chain, not reduced)
                            timing_results['original_lxmixer'] = variant_time
                
                # Use the "optimal" case (not chain, not reduced)
                results['original_lxmixer'] = cnots
                print(f"Original - optimal cost: {cnots}, time: {timing_results['original_lxmixer']:.4f}s")
                
            except Exception as e:
                print(f"Error with Original LXMixer: {e}")
                import traceback
                traceback.print_exc()
                results['original_lxmixer'] = float('inf')
                timing_results['original_lxmixer'] = float('inf')
        else:
            # Skip Original LXMixer - use default value
            results['original_lxmixer'] = 0  # or float('inf') if you prefer
            timing_results['original_lxmixer'] = 0.0
        
        # LXMixer (uses integers) - run for each method
        for method in self.lxmixer_methods:
            method_key = f'lxmixer_{method}'
            timing_key = f'lxmixer_{method}_time'
            
            try:
                start_time = time.time()
                mixer = LXMixer(B_integers, nL, method=method)
                mixer.compute_family_of_valid_graphs()
                mixer.compute_all_orbits()
                mixer.compute_minimal_generating_sets()
                mixer.compute_projector_stabilizers()
                mixer.compute_costs()
                
                mixer.find_best_mixer()
                best_cost = mixer.best_cost
                
                end_time = time.time()
                
                results[method_key] = best_cost
                timing_results[timing_key] = end_time - start_time
                print(f"{method.replace("_", " ").capitalize()} - optimal cost: {best_cost}, time: {timing_results[timing_key]:.4f}s")
                
            except Exception as e:
                print(f"Error with {method.replace("_", " ").capitalize()}: {e}")
                failed_B_info = {
                    'B_strings': B_strings,
                    'B_integers': B_integers,
                    'method': method,
                    'error': str(e)
                }
                results[method_key] = float('inf')
                timing_results[timing_key] = float('inf')
                # Return the failed B info as additional data
                return (*[results.get(key, 0) for key in ['original_lxmixer'] + [f'lxmixer_{m}' for m in self.lxmixer_methods]], 
                        *[timing_results.get(key, 0) for key in ['original_lxmixer'] + [f'lxmixer_{m}_time' for m in self.lxmixer_methods]], 
                        failed_B_info)

        # Print summary for all methods
        summary_parts = [f"Original: {results['original_lxmixer']}"]
        for method in self.lxmixer_methods:
            summary_parts.append(f"{method.replace("_", " ").capitalize()}: {results[f'lxmixer_{method}']}")
        print(", ".join(summary_parts))
        
        return (*[results.get(key, 0) for key in ['original_lxmixer'] + [f'lxmixer_{m}' for m in self.lxmixer_methods]], 
                *[timing_results.get(key, 0) for key in ['original_lxmixer'] + [f'lxmixer_{m}_time' for m in self.lxmixer_methods]], 
                None)


def saveResult(res):
    global cnots_original_lxmixer, cnots_lxmixer_methods, times_original_lxmixer, times_lxmixer_methods, failed_B_sets, infinity_B_sets
    
    # Get the active LXMixer methods for this run
    lxmixer_methods = []
    if RUN_LXMIXER_LARGEST_ORBIT:
        lxmixer_methods.append("largest_orbits")
    if RUN_LXMIXER_ALL_SUBORBIT:
        lxmixer_methods.append("all_suborbits")
    if RUN_LXMIXER_SEMI_RESTRICTED_SUBORBIT:
        lxmixer_methods.append("semi_restricted_suborbits")
    
    num_methods = len(lxmixer_methods)
    
    # Extract Original LXMixer cost (always at index 0)
    cnots_original_lxmixer.append(res[0])
    
    # Extract LXMixer costs for each method (indices 1 to num_methods)
    for i, method in enumerate(lxmixer_methods, 1):
        if method not in cnots_lxmixer_methods:
            cnots_lxmixer_methods[method] = []
        
        if i < len(res) and isinstance(res[i], (int, float)):
            cnots_lxmixer_methods[method].append(res[i])
        else:
            cnots_lxmixer_methods[method].append(float('inf'))
    
    # Extract timing data (starts after costs)
    original_lxmixer_time_idx = 1 + num_methods
    times_original_lxmixer.append(res[original_lxmixer_time_idx] if original_lxmixer_time_idx < len(res) else 0.0)
    
    # Extract LXMixer timing for each method
    for i, method in enumerate(lxmixer_methods):
        if method not in times_lxmixer_methods:
            times_lxmixer_methods[method] = []
        
        time_idx = original_lxmixer_time_idx + 1 + i
        if time_idx < len(res) and isinstance(res[time_idx], (int, float)):
            times_lxmixer_methods[method].append(res[time_idx])
        else:
            times_lxmixer_methods[method].append(0.0)
    
    # Check for failed B set (last element)
    failed_info_idx = original_lxmixer_time_idx + 1 + num_methods
    if failed_info_idx < len(res) and res[failed_info_idx] is not None:
        failed_B_sets.append(res[failed_info_idx])


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
    
    # Build list of active LXMixer methods
    lxmixer_methods = []
    if RUN_LXMIXER_LARGEST_ORBIT:
        lxmixer_methods.append("largest_orbits")
    if RUN_LXMIXER_ALL_SUBORBIT:
        lxmixer_methods.append("all_suborbits")
    if RUN_LXMIXER_SEMI_RESTRICTED_SUBORBIT:
        lxmixer_methods.append("semi_restricted_suborbits")
    
    # Configuration summary
    config_parts = [f"Original={'ON' if RUN_ORIGINAL_LXMIXER else 'OFF'}"]
    for method in lxmixer_methods:
        config_parts.append(f"{method}=ON")
    if not lxmixer_methods:
        config_parts.append("LXMixer=OFF")
    print(f"Configuration: {', '.join(config_parts)}")
    
    # Validate configuration
    if not RUN_ORIGINAL_LXMIXER and not lxmixer_methods:
        print("ERROR: Both mixers are disabled! Please enable at least one mixer.")
        return
    
    worker = Worker(n, checkmixer=False, run_original_lxmixer=RUN_ORIGINAL_LXMIXER, lxmixer_methods=lxmixer_methods)

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

    # Arrays for Original LXMixer
    min_cnots_original_lxmixer = np.zeros(len(m_list))
    max_cnots_original_lxmixer = np.zeros(len(m_list))
    mean_cnots_original_lxmixer = np.zeros(len(m_list))
    var_cnots_original_lxmixer = np.zeros(len(m_list))

    # Arrays for each LXMixer method
    min_cnots_lxmixer_methods = {}
    max_cnots_lxmixer_methods = {}
    mean_cnots_lxmixer_methods = {}
    var_cnots_lxmixer_methods = {}
    
    for method in lxmixer_methods:
        min_cnots_lxmixer_methods[method] = np.zeros(len(m_list))
        max_cnots_lxmixer_methods[method] = np.zeros(len(m_list))
        mean_cnots_lxmixer_methods[method] = np.zeros(len(m_list))
        var_cnots_lxmixer_methods[method] = np.zeros(len(m_list))

    # Arrays for timing data
    min_times_original_lxmixer = np.zeros(len(m_list))
    max_times_original_lxmixer = np.zeros(len(m_list))
    mean_times_original_lxmixer = np.zeros(len(m_list))
    var_times_original_lxmixer = np.zeros(len(m_list))

    min_times_lxmixer_methods = {}
    max_times_lxmixer_methods = {}
    mean_times_lxmixer_methods = {}
    var_times_lxmixer_methods = {}
    
    for method in lxmixer_methods:
        min_times_lxmixer_methods[method] = np.zeros(len(m_list))
        max_times_lxmixer_methods[method] = np.zeros(len(m_list))
        mean_times_lxmixer_methods[method] = np.zeros(len(m_list))
        var_times_lxmixer_methods[method] = np.zeros(len(m_list))

    # Track failed B sets for LXMixer
    all_failed_B_sets = []
    all_infinity_B_sets = []

    for i, m in enumerate(m_list):
        print(f"\nProcessing m={m} (|B|={m}) - {i+1}/{len(m_list)}")
        global cnots_original_lxmixer, cnots_lxmixer_methods, times_original_lxmixer, times_lxmixer_methods, failed_B_sets, infinity_B_sets
        cnots_original_lxmixer = []
        cnots_lxmixer_methods = {method: [] for method in lxmixer_methods}
        times_original_lxmixer = []
        times_lxmixer_methods = {method: [] for method in lxmixer_methods}
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
        
        # Show collected results
        lxmixer_result_counts = {method: len(results) for method, results in cnots_lxmixer_methods.items()}
        print(f"Collected {len(cnots_original_lxmixer)} Original results and {lxmixer_result_counts} LXMixer results for m={m}")
        print(f"Failed: {failed_count}/{num_samples}")
        
        # Debug: Show what values we actually collected
        print(f"Original results: {cnots_original_lxmixer}")
        for method, results in cnots_lxmixer_methods.items():
            print(f"{method.replace("_", " ").capitalize()} results: {results}")
        
        # Store failed B sets for this m value
        if failed_B_sets:
            print(f"LXMixer failed on {len(failed_B_sets)} B sets for m={m}")
            all_failed_B_sets.extend([(m, failed_B) for failed_B in failed_B_sets])
        
        # Check for infinity results and store those B sets
        for method, results in cnots_lxmixer_methods.items():
            for j, lxmixer_result in enumerate(results):
                if lxmixer_result == float('inf'):
                    # Find the corresponding B set
                    if j < len(worker_B_sets):
                        B_strings, B_integers = worker_B_sets[j]
                        infinity_B_info = {
                            'B_strings': B_strings,
                            'B_integers': B_integers,
                            'method': method,
                            'error': f'LXMixer {method.replace("_", " ").capitalize()} returned infinity (no exception thrown)'
                        }
                        infinity_B_sets.append(infinity_B_info)
        
        if infinity_B_sets:
            print(f"LXMixer returned infinity for {len(infinity_B_sets)} B sets for m={m}")
            all_infinity_B_sets.extend([(m, infinity_B) for infinity_B in infinity_B_sets])
        
        # Check if we have any valid data
        if len(cnots_original_lxmixer) == 0 or all(len(results) == 0 for results in cnots_lxmixer_methods.values()):
            print(f"No results collected for m={m} - all tasks failed!")
            # Set to NaN instead of skipping to maintain array indices
            min_cnots_original_lxmixer[i] = np.nan
            max_cnots_original_lxmixer[i] = np.nan
            mean_cnots_original_lxmixer[i] = np.nan
            var_cnots_original_lxmixer[i] = np.nan
            
            for method in lxmixer_methods:
                min_cnots_lxmixer_methods[method][i] = np.nan
                max_cnots_lxmixer_methods[method][i] = np.nan
                mean_cnots_lxmixer_methods[method][i] = np.nan
                var_cnots_lxmixer_methods[method][i] = np.nan
                
                min_times_lxmixer_methods[method][i] = np.nan
                max_times_lxmixer_methods[method][i] = np.nan
                mean_times_lxmixer_methods[method][i] = np.nan
                var_times_lxmixer_methods[method][i] = np.nan
            
            # Timing data
            min_times_original_lxmixer[i] = np.nan
            max_times_original_lxmixer[i] = np.nan
            mean_times_original_lxmixer[i] = np.nan
            var_times_original_lxmixer[i] = np.nan
            continue

        # Calculate statistics for Original LXMixer
        min_cnots_original_lxmixer[i] = np.min(cnots_original_lxmixer)
        max_cnots_original_lxmixer[i] = np.max(cnots_original_lxmixer)
        mean_cnots_original_lxmixer[i] = np.mean(cnots_original_lxmixer)
        var_cnots_original_lxmixer[i] = np.var(cnots_original_lxmixer)

        # Calculate statistics for each LXMixer method
        for method in lxmixer_methods:
            method_results = cnots_lxmixer_methods[method]
            print(f"{method.replace("_", " ").capitalize()} results for m={m}: {method_results}")
            finite_results = [x for x in method_results if np.isfinite(x)]
            print(f"{method.replace("_", " ").capitalize()} finite values: {finite_results}")
            
            if len(finite_results) == 0:
                print(f"All {method.replace("_", " ").capitalize()} values are infinite/NaN for m={m}")
                min_cnots_lxmixer_methods[method][i] = np.nan
                max_cnots_lxmixer_methods[method][i] = np.nan
                mean_cnots_lxmixer_methods[method][i] = np.nan
                var_cnots_lxmixer_methods[method][i] = np.nan
            else:
                min_cnots_lxmixer_methods[method][i] = np.min(finite_results)
                max_cnots_lxmixer_methods[method][i] = np.max(finite_results)
                mean_cnots_lxmixer_methods[method][i] = np.mean(finite_results)
                var_cnots_lxmixer_methods[method][i] = np.var(finite_results)

        # Calculate timing statistics for Original LXMixer
        min_times_original_lxmixer[i] = np.min(times_original_lxmixer)
        max_times_original_lxmixer[i] = np.max(times_original_lxmixer)
        mean_times_original_lxmixer[i] = np.mean(times_original_lxmixer)
        var_times_original_lxmixer[i] = np.var(times_original_lxmixer)

        # Calculate timing statistics for each LXMixer method
        for method in lxmixer_methods:
            method_times = times_lxmixer_methods[method]
            min_times_lxmixer_methods[method][i] = np.min(method_times)
            max_times_lxmixer_methods[method][i] = np.max(method_times)
            mean_times_lxmixer_methods[method][i] = np.mean(method_times)
            var_times_lxmixer_methods[method][i] = np.var(method_times)

        print(int(100 * (i + 1) / len(m_list)), "%")

    # Create backward compatibility variables for plotting (combine all LXMixer methods)
    # This is a temporary solution until we update the plotting code
    if lxmixer_methods:
        # For plotting, we'll use the first active method as representative
        first_method = lxmixer_methods[0]
        mean_cnots_lxmixer = mean_cnots_lxmixer_methods[first_method].copy()
        var_cnots_lxmixer = var_cnots_lxmixer_methods[first_method].copy()
        min_cnots_lxmixer = min_cnots_lxmixer_methods[first_method].copy()
        max_cnots_lxmixer = max_cnots_lxmixer_methods[first_method].copy()
        mean_times_lxmixer = mean_times_lxmixer_methods[first_method].copy()
        var_times_lxmixer = var_times_lxmixer_methods[first_method].copy()
        min_times_lxmixer = min_times_lxmixer_methods[first_method].copy()
        max_times_lxmixer = max_times_lxmixer_methods[first_method].copy()
    else:
        # No LXMixer methods enabled
        mean_cnots_lxmixer = np.full(len(m_list), np.nan)
        var_cnots_lxmixer = np.full(len(m_list), np.nan)
        min_cnots_lxmixer = np.full(len(m_list), np.nan)
        max_cnots_lxmixer = np.full(len(m_list), np.nan)
        mean_times_lxmixer = np.full(len(m_list), np.nan)
        var_times_lxmixer = np.full(len(m_list), np.nan)
        min_times_lxmixer = np.full(len(m_list), np.nan)
        max_times_lxmixer = np.full(len(m_list), np.nan)

    # Print final data summary
    print(f"\nFinal data summary:")
    print(f"m_list: {m_list}")
    
    if RUN_ORIGINAL_LXMIXER:
        print(f"Original means: {mean_cnots_original_lxmixer}")
        print(f"Original time means: {mean_times_original_lxmixer}")
        print(f"Original has {np.sum(~np.isnan(mean_cnots_original_lxmixer))} valid values")
    else:
        print("Original was disabled")
    
    if lxmixer_methods:
        for method in lxmixer_methods:
            print(f"{method.replace("_", " ").capitalize()} means: {mean_cnots_lxmixer_methods[method]}")
            print(f"{method.replace("_", " ").capitalize()} time means: {mean_times_lxmixer_methods[method]}")
            print(f"{method.replace("_", " ").capitalize()} has {np.sum(~np.isnan(mean_cnots_lxmixer_methods[method]))} valid values")
    else:
        print("LXMixer was disabled")

    # Debug: Check for infinite values
    print(f"\nDebugging infinite values:")
    if RUN_ORIGINAL_LXMIXER:
        print(f"Original infinite values: {np.sum(np.isinf(mean_cnots_original_lxmixer))}")
    if lxmixer_methods:
        for method in lxmixer_methods:
            print(f"{method.replace("_", " ").capitalize()} infinite values: {np.sum(np.isinf(mean_cnots_lxmixer_methods[method]))}")
    
    # Check each m value individually
    for i, m in enumerate(m_list):
        debug_parts = []
        
        if RUN_ORIGINAL_LXMIXER:
            original_lxmixer_val = mean_cnots_original_lxmixer[i]
            debug_parts.append(f"Original={original_lxmixer_val}")
        
        # Add each LXMixer method individually
        for method in lxmixer_methods:
            method_val = mean_cnots_lxmixer_methods[method][i]
            debug_parts.append(f"{method.replace("_", " ").capitalize()}={method_val}")
        
        if debug_parts:
            print(f"m={m}: {', '.join(debug_parts)}")
            
            # Check for NaN/Infinity in all values
            all_values = []
            if RUN_ORIGINAL_LXMIXER:
                all_values.append(mean_cnots_original_lxmixer[i])
            for method in lxmixer_methods:
                all_values.append(mean_cnots_lxmixer_methods[method][i])
            
            if any(np.isnan(val) for val in all_values):
                print(f"  -> NaN detected")
            if any(np.isinf(val) for val in all_values):
                print(f"  -> Infinity detected")

    # Create comparison plot for CNOT costs and timing
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Filter out NaN and infinite values for plotting based on enabled mixers
    if RUN_ORIGINAL_LXMIXER and RUN_LXMIXER:
        # Both mixers enabled - filter based on both
        valid_indices = ~(np.isnan(mean_cnots_original_lxmixer) | np.isnan(mean_cnots_lxmixer) | 
                          np.isinf(mean_cnots_original_lxmixer) | np.isinf(mean_cnots_lxmixer))
    elif RUN_ORIGINAL_LXMIXER:
        # Only Original LXMixer enabled
        valid_indices = ~(np.isnan(mean_cnots_original_lxmixer) | np.isinf(mean_cnots_original_lxmixer))
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
    if RUN_ORIGINAL_LXMIXER:
        mean_original_lxmixer_valid = mean_cnots_original_lxmixer[valid_indices]
        var_original_lxmixer_valid = var_cnots_original_lxmixer[valid_indices]
        min_original_lxmixer_valid = min_cnots_original_lxmixer[valid_indices]
        max_original_lxmixer_valid = max_cnots_original_lxmixer[valid_indices]
    
    if RUN_LXMIXER:
        mean_lxmixer_valid = mean_cnots_lxmixer[valid_indices]
        var_lxmixer_valid = var_cnots_lxmixer[valid_indices]
        min_lxmixer_valid = min_cnots_lxmixer[valid_indices]
        max_lxmixer_valid = max_cnots_lxmixer[valid_indices]

    # Timing data
    if RUN_ORIGINAL_LXMIXER:
        mean_times_original_lxmixer_valid = mean_times_original_lxmixer[valid_indices]
        var_times_original_lxmixer_valid = var_times_original_lxmixer[valid_indices]
        min_times_original_lxmixer_valid = min_times_original_lxmixer[valid_indices]
        max_times_original_lxmixer_valid = max_times_original_lxmixer[valid_indices]
    
    if RUN_LXMIXER:
        mean_times_lxmixer_valid = mean_times_lxmixer[valid_indices]
        var_times_lxmixer_valid = var_times_lxmixer[valid_indices]
        min_times_lxmixer_valid = min_times_lxmixer[valid_indices]
        max_times_lxmixer_valid = max_times_lxmixer[valid_indices]

    # Plot CNOT costs (left subplot)
    if RUN_ORIGINAL_LXMIXER:
        ax1.plot(m_list_valid, mean_original_lxmixer_valid, "o-b", label="Original (mean)")
        ax1.fill_between(
            m_list_valid,
            mean_original_lxmixer_valid - np.sqrt(var_original_lxmixer_valid),
            mean_original_lxmixer_valid + np.sqrt(var_original_lxmixer_valid),
            color="blue",
            alpha=0.35,
            label="Original (1σ)"
        )
        ax1.plot(m_list_valid, min_original_lxmixer_valid, ":b", alpha=0.7)
        ax1.plot(m_list_valid, max_original_lxmixer_valid, ":b", alpha=0.7)

    # Plot each LXMixer method separately
    colors = ["red", "green", "orange", "purple", "brown"]
    markers = ["x", "s", "^", "v", "d"]
    
    for i, method in enumerate(lxmixer_methods):
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        
        method_mean_valid = mean_cnots_lxmixer_methods[method][valid_indices]
        method_var_valid = var_cnots_lxmixer_methods[method][valid_indices]
        method_min_valid = min_cnots_lxmixer_methods[method][valid_indices]
        method_max_valid = max_cnots_lxmixer_methods[method][valid_indices]
        
        ax1.plot(m_list_valid, method_mean_valid, f"{marker}-", color=color, label=f"{method.replace("_", " ").capitalize()} (mean)")
        ax1.fill_between(
            m_list_valid,
            method_mean_valid - np.sqrt(method_var_valid),
            method_mean_valid + np.sqrt(method_var_valid),
            color=color,
            alpha=0.35,
            label=f"{method.replace("_", " ").capitalize()} (1σ)"
        )
        ax1.plot(m_list_valid, method_min_valid, ":", color=color, alpha=0.7)
        ax1.plot(m_list_valid, method_max_valid, ":", color=color, alpha=0.7)

    ax1.legend()
    ax1.set_xlim(min(m_list_valid), max(m_list_valid))
    ax1.set_xlabel(r"$|B|$")
    ax1.set_ylabel("#CNOTS")
    
    # Update title based on enabled mixers
    if RUN_ORIGINAL_LXMIXER and lxmixer_methods:
        ax1.set_title(f"CNOT Cost Comparison (n={n})")
    elif RUN_ORIGINAL_LXMIXER:
        ax1.set_title(f"Original CNOT Cost (n={n})")
    else:
        method_names = ", ".join(lxmixer_methods)
        ax1.set_title(f"{method_names.replace("_", " ").capitalize()} CNOT Cost () (n={n})")
    ax1.grid()

    # Plot execution times (right subplot)
    if RUN_ORIGINAL_LXMIXER:
        ax2.plot(m_list_valid, mean_times_original_lxmixer_valid, "o-b", label="Original (mean)")
        ax2.fill_between(
            m_list_valid,
            mean_times_original_lxmixer_valid - np.sqrt(var_times_original_lxmixer_valid),
            mean_times_original_lxmixer_valid + np.sqrt(var_times_original_lxmixer_valid),
            color="blue",
            alpha=0.35,
            label="Original (1σ)"
        )
        ax2.plot(m_list_valid, min_times_original_lxmixer_valid, ":b", alpha=0.7)
        ax2.plot(m_list_valid, max_times_original_lxmixer_valid, ":b", alpha=0.7)

    # Plot timing for each LXMixer method separately
    for i, method in enumerate(lxmixer_methods):
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        
        method_mean_times_valid = mean_times_lxmixer_methods[method][valid_indices]
        method_var_times_valid = var_times_lxmixer_methods[method][valid_indices]
        method_min_times_valid = min_times_lxmixer_methods[method][valid_indices]
        method_max_times_valid = max_times_lxmixer_methods[method][valid_indices]
        
        ax2.plot(m_list_valid, method_mean_times_valid, f"{marker}-", color=color, label=f"{method.replace("_", " ").capitalize()} (mean)")
        ax2.fill_between(
            m_list_valid,
            method_mean_times_valid - np.sqrt(method_var_times_valid),
            method_mean_times_valid + np.sqrt(method_var_times_valid),
            color=color,
            alpha=0.35,
            label=f"{method.replace("_", " ").capitalize()} (1σ)"
        )
        ax2.plot(m_list_valid, method_min_times_valid, ":", color=color, alpha=0.7)
        ax2.plot(m_list_valid, method_max_times_valid, ":", color=color, alpha=0.7)

    ax2.legend()
    ax2.set_xlim(min(m_list_valid), max(m_list_valid))
    ax2.set_xlabel(r"$|B|$")
    ax2.set_ylabel("Execution Time (s)")
    
    # Update title based on enabled mixers
    if RUN_ORIGINAL_LXMIXER and lxmixer_methods:
        ax2.set_title(f"Execution Time Comparison (n={n})")
    elif RUN_ORIGINAL_LXMIXER:
        ax2.set_title(f"Original Execution Time (n={n})")
    else:
        method_names = ", ".join(lxmixer_methods)
        ax2.set_title(f"{method_names.replace("_", " ").capitalize()} Execution Time (n={n})")
    ax2.grid()

    plt.tight_layout()
    
    # Update filename based on enabled mixers
    if RUN_ORIGINAL_LXMIXER and lxmixer_methods:
        plt.savefig(f"comparison_n{n}_multi_methods.png")
    elif RUN_ORIGINAL_LXMIXER:
        plt.savefig(f"Original_LXMixer_only_n{n}.png")
    else:
        method_str = "_".join(lxmixer_methods)
        plt.savefig(f"LXMixer_{method_str}_n{n}.png")
    plt.clf()

    # Create separate timing plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if RUN_ORIGINAL_LXMIXER:
        ax.plot(m_list_valid, mean_times_original_lxmixer_valid, "o-b", label="Original (mean)")
        ax.fill_between(
            m_list_valid,
            mean_times_original_lxmixer_valid - np.sqrt(var_times_original_lxmixer_valid),
            mean_times_original_lxmixer_valid + np.sqrt(var_times_original_lxmixer_valid),
            color="blue",
            alpha=0.35,
            label="Original (1σ)"
        )
        ax.plot(m_list_valid, min_times_original_lxmixer_valid, ":b", alpha=0.7)
        ax.plot(m_list_valid, max_times_original_lxmixer_valid, ":b", alpha=0.7)

    # Plot timing for each LXMixer method separately
    colors = ["red", "green", "orange", "purple", "brown"]
    markers = ["x", "s", "^", "v", "d"]
    
    for i, method in enumerate(lxmixer_methods):
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        
        method_mean_times_valid = mean_times_lxmixer_methods[method][valid_indices]
        method_var_times_valid = var_times_lxmixer_methods[method][valid_indices]
        method_min_times_valid = min_times_lxmixer_methods[method][valid_indices]
        method_max_times_valid = max_times_lxmixer_methods[method][valid_indices]
        
        ax.plot(m_list_valid, method_mean_times_valid, f"{marker}-", color=color, label=f"{method.replace("_", " ").capitalize()} (mean)")
        ax.fill_between(
            m_list_valid,
            method_mean_times_valid - np.sqrt(method_var_times_valid),
            method_mean_times_valid + np.sqrt(method_var_times_valid),
            color=color,
            alpha=0.35,
            label=f"{method.replace("_", " ").capitalize()} (1σ)"
        )
        ax.plot(m_list_valid, method_min_times_valid, ":", color=color, alpha=0.7)
        ax.plot(m_list_valid, method_max_times_valid, ":", color=color, alpha=0.7)

    ax.legend()
    ax.set_xlim(min(m_list_valid), max(m_list_valid))
    ax.set_xlabel(r"$|B|$")
    ax.set_ylabel("Execution Time (s)")
    
    # Update title based on enabled mixers
    if RUN_ORIGINAL_LXMIXER and lxmixer_methods:
        ax.set_title(f"Execution Time Comparison: Original LXMixer vs Individual LXMixer Methods (n={n})")
        timing_filename = f"timing_comparison_n{n}.png"
    elif RUN_ORIGINAL_LXMIXER:
        ax.set_title(f"Original LXMixer Execution Time (n={n})")
        timing_filename = f"timing_original_lxmixer_only_n{n}.png"
    else:
        method_names = ", ".join([method.replace("_", " ").capitalize() for method in lxmixer_methods])
        ax.set_title(f"LXMixer Execution Time: {method_names} (n={n})")
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
    
    # How to use the configuration:
    # 
    # Run selected mixers (default configuration above, set 
    # RUN_ORIGINAL_LXMIXER = True/False
    # RUN_LXMIXER_LARGEST_ORBIT = True/False       
    # RUN_LXMIXER_ALL_SUBORBIT = True/False         
    # RUN_LXMIXER_SEMI_RESTRICTED_SUBORBIT = False (not implemented yet)
    #
    #    for n in [3]:
    #        main(n)
    
    for n in [3]:
        main(n)