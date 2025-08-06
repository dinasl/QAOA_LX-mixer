import time
import numpy as np
from random import sample
import itertools
import sys
import os

# Add parent directory to path to import Mixer
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mixer_path = os.path.join(parent_dir, 'Mixer.py')
sys.path.append(parent_dir)

from Mixer import *
from Stabilizer import *

# Importing Franz mixer and related modules
sys.path.append(r"C:\Users\sanne\LogicalXMixer")
from Mixer_Franz import MixerFranz
from PauliString import *
from PauliOperations import *
from Stabilizers import *
from GroupGraph import *
import networkx as nx

def generate_random_B(n, m):
    """
    Generate a random B set of size m for n qubits.
    
    Args:
        n (int): Number of qubits
        m (int): Size of B set (|B|)
        
    Returns:
        tuple: (B_strings, B_integers) - B set in both string and integer formats
    """
    all_states = ["".join(i) for i in itertools.product("01", repeat=n)]
    binary_strings = sample(all_states, m)
    binary_integers = [int(binary_str, 2) for binary_str in binary_strings]
    return binary_strings, binary_integers

class MixerAnalyzer:
    """
    Enhanced analyzer for Franz mixer with detailed internal analysis.
    """
    def __init__(self, B, reduced=True, digraph=True, sort=True):
        """
        Initialize the mixer analyzer.
        
        Args:
            B: List of binary strings representing the basis states
            reduced: Whether to use reduced projectors
            digraph: Whether to use directed graphs
            sort: Whether to sort the basis states
        """
        self.B = B
        self.reduced = reduced
        self.digraph = digraph
        self.sort = sort
        
        # Initialize the mixer
        self.mixer = MixerFranz(B, digraph=digraph, reduced=reduced, sort=sort)
        
        # Find the best mixer
        self.mixer.get_best_mixer_commuting_graphs()
        
        # Store the solution
        if reduced:
            self.solution = self.mixer.solution_reduced[0] if hasattr(self.mixer, 'solution_reduced') else []
            self.total_cost = self.mixer.solution_reduced_cost if hasattr(self.mixer, 'solution_reduced_cost') else 0
        else:
            self.solution = self.mixer.solution[0] if hasattr(self.mixer, 'solution') else []
            self.total_cost = self.mixer.solution_cost if hasattr(self.mixer, 'solution_cost') else 0

    def analyze_solution(self):
        """
        Perform comprehensive analysis of the optimal solution.
        """
        print("  DETAILED FRANZ MIXER ANALYSIS:")
        print("  " + "-"*50)
        
        print(f"  Basis states B: {[b.state if hasattr(b, 'state') else str(b) for b in self.mixer.B]}")
        print(f"  Number of basis states: {len(self.mixer.B)}")
        print(f"  Qubit count: {len(self.mixer.B[0].state) if hasattr(self.mixer.B[0], 'state') else len(str(self.mixer.B[0]))}")
        print(f"  Using reduced projectors: {self.reduced}")
        print(f"  Total optimal cost: {self.total_cost}")
        print(f"  Number of graphs in solution: {len(self.solution)}")
        print()
        
        if not self.solution:
            print("  No solution found!")
            return
        
        # Analyze each graph in the solution
        for i, graph in enumerate(self.solution):
            self._analyze_single_graph(i, graph)
        
        # Overall connectivity analysis
        self._analyze_connectivity()
        
        # Cost breakdown
        self._analyze_cost_breakdown()

    def _analyze_single_graph(self, index, graph):
        """
        Analyze a single graph in the solution.
        """
        print(f"  GRAPH {index + 1}:")
        print("  " + "-" * 30)
        
        # Basic graph info
        if hasattr(graph, 'Xl'):
            print(f"    Logical X operator: {graph.Xl}")
        
        if self.reduced:
            cost = getattr(graph, 'cost_reduced', 'N/A')
            ps_count = len(getattr(graph, 'PS_reduced', []))
            print(f"    Cost (reduced): {cost}")
            print(f"    Pauli strings (reduced): {ps_count} terms")
        else:
            cost = getattr(graph, 'cost', 'N/A')
            ps_count = len(getattr(graph, 'PS', []))
            print(f"    Cost (full): {cost}")
            print(f"    Pauli strings (full): {ps_count} terms")
        
        # Graph structure
        if hasattr(graph, 'G'):
            print(f"    Graph edges: {list(graph.G.edges)}")
            print(f"    Graph nodes: {list(graph.G.nodes)}")
            print(f"    Number of edges: {graph.G.number_of_edges()}")
            print(f"    Number of nodes: {graph.G.number_of_nodes()}")
        
        # Show the projector terms (limited to first few for readability)
        if self.reduced and hasattr(graph, 'PS_reduced'):
            print("    Reduced Projector Terms (first 5):")
            for j, ps in enumerate(graph.PS_reduced[:5]):
                if hasattr(ps, 'c') and hasattr(ps, 'P'):
                    print(f"      Term {j+1}: {ps.c:.4f} * {ps.P}")
                else:
                    print(f"      Term {j+1}: {ps}")
            if len(graph.PS_reduced) > 5:
                print(f"      ... and {len(graph.PS_reduced) - 5} more terms")
        elif not self.reduced and hasattr(graph, 'PS'):
            print("    Full Projector Terms (first 5):")
            for j, ps in enumerate(graph.PS[:5]):
                if hasattr(ps, 'c') and hasattr(ps, 'P'):
                    print(f"      Term {j+1}: {ps.c:.4f} * {ps.P}")
                else:
                    print(f"      Term {j+1}: {ps}")
            if len(graph.PS) > 5:
                print(f"      ... and {len(graph.PS) - 5} more terms")
        
        print()

    def _analyze_connectivity(self):
        """
        Analyze the connectivity of the combined graph.
        """
        print("  CONNECTIVITY ANALYSIS:")
        print("  " + "-" * 30)
        
        try:
            # Combine all graphs
            if self.digraph:
                combined_G = nx.DiGraph()
            else:
                combined_G = nx.Graph()
            
            for graph in self.solution:
                if hasattr(graph, 'G'):
                    combined_G = nx.compose(combined_G, graph.G)
            
            print(f"    Combined graph nodes: {list(combined_G.nodes)}")
            print(f"    Combined graph edges: {list(combined_G.edges)}")
            print(f"    Total edges: {combined_G.number_of_edges()}")
            
            # Check connectivity
            if self.digraph:
                is_connected = nx.is_connected(combined_G.to_undirected())
            else:
                is_connected = nx.is_connected(combined_G)
            
            print(f"    Is connected: {is_connected}")
            
            if is_connected and len(combined_G.nodes) > 1:
                if self.digraph:
                    diameter = nx.diameter(combined_G.to_undirected())
                else:
                    diameter = nx.diameter(combined_G)
                print(f"    Graph diameter: {diameter}")
        
        except Exception as e:
            print(f"    Error in connectivity analysis: {e}")
        
        print()

    def _analyze_cost_breakdown(self):
        """
        Analyze the cost breakdown.
        """
        print("  COST BREAKDOWN:")
        print("  " + "-" * 30)
        
        individual_costs = []
        for i, graph in enumerate(self.solution):
            if self.reduced:
                cost = getattr(graph, 'cost_reduced', 0)
            else:
                cost = getattr(graph, 'cost', 0)
            individual_costs.append(cost)
            print(f"    Graph {i+1} cost: {cost}")
        
        print(f"    Sum of individual costs: {sum(individual_costs)}")
        print(f"    Total reported cost: {self.total_cost}")
        
        if hasattr(self.mixer, 'base_cost'):
            print(f"    Base cost: {self.mixer.base_cost}")
        print()

    def get_stabilizer_info(self):
        """
        Extract stabilizer generator information where possible.
        """
        try:
            stabilizer_info = []
            
            for i, graph in enumerate(self.solution):
                graph_info = {
                    'graph_index': i,
                    'logical_X': getattr(graph, 'Xl', 'N/A'),
                    'cost': getattr(graph, 'cost_reduced' if self.reduced else 'cost', 'N/A'),
                    'stabilizer_generators': 'Analysis requires additional Franz mixer modules'
                }
                stabilizer_info.append(graph_info)
            
            return stabilizer_info
            
        except Exception as e:
            return [{'error': f'Error extracting stabilizer info: {e}'}]

    def compare_with_chain_mixer(self):
        """
        Compare the optimal solution with a chain mixer.
        """
        print("  COMPARISON WITH CHAIN MIXER:")
        print("  " + "-" * 30)
        
        try:
            # Get chain mixer results
            self.mixer.get_chain_mixer()
            
            if self.reduced:
                chain_cost = getattr(self.mixer, 'solution_chain_reduced_cost', float('inf'))
                chain_solution = getattr(self.mixer, 'solution_chain_reduced', [])
            else:
                chain_cost = getattr(self.mixer, 'solution_chain_cost', float('inf'))
                chain_solution = getattr(self.mixer, 'solution_chain', [])
            
            print(f"    Optimal mixer cost: {self.total_cost}")
            print(f"    Chain mixer cost: {chain_cost}")
            print(f"    Improvement: {chain_cost - self.total_cost} gates")
            if self.total_cost > 0:
                print(f"    Improvement ratio: {chain_cost / self.total_cost:.2f}x")
            print(f"    Chain mixer graphs: {len(chain_solution)}")
            
        except Exception as e:
            print(f"    Error in chain mixer comparison: {e}")
        
        print()

def analyze_lxmixer(B_integers, n):
    """
    Analyze LXMixer in detail for a given B set.
    
    Args:
        B_integers (list): B set as integers
        n (int): Number of qubits
        
    Returns:
        dict: Detailed analysis results
    """
    print("="*60)
    print("LXMixer Analysis")
    print("="*60)
    
    results = {}
    
    try:
        # Initialize LXMixer
        start_time = time.time()
        mixer = LXMixer(B_integers, n)
        
        # Compute family of valid graphs
        print("Computing family of valid graphs...")
        mixer.compute_family_of_valid_graphs()
        
        print(f"Family of valid graphs ({len(mixer.family_of_valid_graphs)} operators):")
        for X, edges in mixer.family_of_valid_graphs.items():
            print(f"  X_{X:0{n}b}: {edges}")
        
        # Compute all orbits
        print("\nComputing all orbits...")
        mixer.compute_all_orbits()
        
        print(f"\nOrbits found ({len(mixer.orbits)} orbits):")
        for i, (nodes, orbit) in enumerate(mixer.orbits.items()):
            print(f"  Orbit {i+1}: Nodes {nodes}")
            print(f"    Logical X operators: [{', '.join(f'{X:0{n}b}' for X in orbit.Xs)}]")
        
        # Initialize stabilizer and compute minimal generating sets
        print("\nComputing minimal generating sets...")
        stabilizer = Stabilizer(B=mixer.B, n=mixer.nL, orbit_dictionary=mixer.orbits)
        stabilizer.compute_minimal_generating_sets()
        
        # Compute projectors
        print("Computing projector stabilizers...")
        stabilizer.compute_projector_stabilizers()
        
        # Compute costs
        print("Computing costs...")
        mixer.compute_costs()
        
        print(f"\nDetailed orbit analysis:")
        for i, (nodes, orbit) in enumerate(mixer.orbits.items()):
            print(f"  Orbit {i+1}: Nodes {nodes}")
            print(f"    Logical X operators: [{', '.join(f'{X:0{n}b}' for X in orbit.Xs)}]")
            print(f"    Projectors (Z ops): [{', '.join(f'{"+" if Z[0] == 1 else "-"}{Z[1]:0{n}b}' for Z in orbit.Zs)}]")
            print(f"    Cost: {orbit.cost}")
        
        # Find best mixer
        print("\nFinding best mixer combination...")
        mixer.find_best_mixer()
        
        end_time = time.time()
        
        # Store results
        results = {
            'success': True,
            'execution_time': end_time - start_time,
            'family_of_valid_graphs': mixer.family_of_valid_graphs,
            'orbits': mixer.orbits,
            'best_combinations': mixer.best_combinations,
            'best_Xs': mixer.best_Xs,
            'best_Zs': mixer.best_Zs,
            'best_cost': mixer.best_cost,
            'mixer_object': mixer
        }
        
        print(f"\nBest mixer found:")
        print(f"  Cost: {mixer.best_cost}")
        print(f"  Number of combinations: {len(mixer.best_combinations)}")
        for i, combination in enumerate(mixer.best_combinations):
            print(f"  Combination {i+1}: {combination}")
        
        print(f"\nExecution time: {results['execution_time']:.4f} seconds")
        
    except Exception as e:
        print(f"Error in LXMixer analysis: {e}")
        import traceback
        traceback.print_exc()
        results = {
            'success': False,
            'error': str(e),
            'execution_time': float('inf'),
            'best_cost': float('inf')
        }
    
    return results

def analyze_franz_mixer(B_strings, n):
    """
    Analyze Franz mixer in detail for a given B set.
    
    Args:
        B_strings (list): B set as binary strings
        n (int): Number of qubits
        
    Returns:
        dict: Detailed analysis results
    """
    print("="*60)
    print("Franz Mixer Analysis")
    print("="*60)
    
    results = {}
    
    try:
        start_time = time.time()
        
        # Test all combinations of chain/non-chain and reduced/non-reduced
        mixer_results = {}
        
        configurations = [
            ('chain_reduced', True, True),
            ('chain_full', True, False),
            ('optimal_reduced', False, True),
            ('optimal_full', False, False)
        ]
        
        # Store the best analyzer for detailed analysis
        best_analyzer = None
        best_cost = float('inf')
        
        for config_name, chain, reduced in configurations:
            print(f"\nTesting {config_name} configuration (chain={chain}, reduced={reduced})...")
            
            try:
                mixer_franz = MixerFranz(B_strings, reduced=reduced)
                
                if chain:
                    mixer_franz.get_chain_mixer()
                    if reduced:
                        cost = mixer_franz.solution_chain_reduced_cost
                    else:
                        cost = mixer_franz.solution_chain_cost
                else:
                    mixer_franz.get_best_mixer_commuting_graphs()
                    if reduced:
                        cost = mixer_franz.solution_reduced_cost
                    else:
                        cost = mixer_franz.solution_cost
                
                mixer_results[config_name] = {
                    'cost': cost,
                    'chain': chain,
                    'reduced': reduced,
                    'mixer_object': mixer_franz
                }
                
                print(f"  Cost: {cost}")
                
                # Create analyzer for the optimal (non-chain) configurations
                if not chain and cost < best_cost:
                    try:
                        best_analyzer = MixerAnalyzer(B_strings, reduced=reduced, digraph=True, sort=True)
                        best_cost = cost
                    except Exception as analyzer_error:
                        print(f"  Warning: Could not create detailed analyzer: {analyzer_error}")
                
                # Try to extract more detailed information if available
                if hasattr(mixer_franz, 'solution_X_operators'):
                    print(f"  X operators: {mixer_franz.solution_X_operators}")
                if hasattr(mixer_franz, 'solution_graphs'):
                    print(f"  Graphs used: {len(mixer_franz.solution_graphs) if mixer_franz.solution_graphs else 0}")
                
            except Exception as config_error:
                print(f"  Error in {config_name} configuration: {config_error}")
                mixer_results[config_name] = {
                    'cost': float('inf'),
                    'error': str(config_error)
                }
        
        # Perform detailed analysis with the best analyzer
        if best_analyzer:
            print(f"\n" + "="*60)
            print("DETAILED ANALYSIS (Best Configuration)")
            print("="*60)
            try:
                best_analyzer.analyze_solution()
                
                # Get stabilizer information
                stabilizer_info = best_analyzer.get_stabilizer_info()
                
                # Compare with chain mixer
                best_analyzer.compare_with_chain_mixer()
                
            except Exception as detail_error:
                print(f"Error in detailed analysis: {detail_error}")
                import traceback
                traceback.print_exc()
        
        end_time = time.time()
        
        # Find the best configuration
        best_config = min(mixer_results.keys(), 
                         key=lambda k: mixer_results[k].get('cost', float('inf')))
        best_cost = mixer_results[best_config].get('cost', float('inf'))
        
        results = {
            'success': True,
            'execution_time': end_time - start_time,
            'configurations': mixer_results,
            'best_configuration': best_config,
            'best_cost': best_cost,
            'detailed_analyzer': best_analyzer,
            'stabilizer_info': best_analyzer.get_stabilizer_info() if best_analyzer else None
        }
        
        print(f"\nBest Franz configuration: {best_config}")
        print(f"Best cost: {best_cost}")
        print(f"Execution time: {results['execution_time']:.4f} seconds")
        
    except Exception as e:
        print(f"Error in Franz mixer analysis: {e}")
        import traceback
        traceback.print_exc()
        results = {
            'success': False,
            'error': str(e),
            'execution_time': float('inf'),
            'best_cost': float('inf')
        }
    
    return results

def compare_mixers(n, m, B_strings=None, B_integers=None):
    """
    Compare LXMixer and Franz mixer for a specific case.
    
    Args:
        n (int): Number of qubits
        m (int): Size of B set (|B|)
        B_strings (list, optional): Predefined B set as strings
        B_integers (list, optional): Predefined B set as integers
    """
    print("="*80)
    print(f"MIXER COMPARISON: n={n}, |B|={m}")
    print("="*80)
    
    # Generate or use provided B set
    if B_strings is None or B_integers is None:
        print(f"Generating random B set with {m} elements...")
        B_strings, B_integers = generate_random_B(n, m)
    else:
        print(f"Using provided B set...")
    
    print(f"B (strings): {B_strings}")
    print(f"B (integers): {B_integers}")
    print(f"B (binary): {[f'{b:0{n}b}' for b in B_integers]}")
    
    # Analyze LXMixer
    print("\n")
    lx_results = analyze_lxmixer(B_integers, n)
    
    # Analyze Franz mixer
    print("\n")
    franz_results = analyze_franz_mixer(B_strings, n)
    
    # Comparison summary
    print("\n" + "="*80)
    print("COMPARISON SUMMARY")
    print("="*80)
    
    print(f"LXMixer:")
    if lx_results['success']:
        print(f"  Status: SUCCESS")
        print(f"  Best cost: {lx_results['best_cost']}")
        print(f"  Execution time: {lx_results['execution_time']:.4f}s")
        print(f"  Number of orbits: {len(lx_results['orbits'])}")
        print(f"  Number of best combinations: {len(lx_results['best_combinations'])}")
        
        # Show detailed orbit information
        print(f"  Orbit details:")
        for i, (nodes, orbit) in enumerate(lx_results['orbits'].items()):
            print(f"    Orbit {i+1}: Nodes {nodes}, Cost {orbit.cost}")
            print(f"      X ops: [{', '.join(f'{X:0{n}b}' for X in orbit.Xs)}]")
            print(f"      Z ops: [{', '.join(f'{"+" if Z[0] == 1 else "-"}{Z[1]:0{n}b}' for Z in orbit.Zs[:3])}{'...' if len(orbit.Zs) > 3 else ''}]")
    else:
        print(f"  Status: FAILED")
        print(f"  Error: {lx_results['error']}")
    
    print(f"\nFranz Mixer:")
    if franz_results['success']:
        print(f"  Status: SUCCESS")
        print(f"  Best cost: {franz_results['best_cost']}")
        print(f"  Best configuration: {franz_results['best_configuration']}")
        print(f"  Execution time: {franz_results['execution_time']:.4f}s")
        
        # Show all configuration results
        print(f"  All configurations:")
        for config, data in franz_results['configurations'].items():
            if 'error' not in data:
                print(f"    {config}: {data['cost']}")
            else:
                print(f"    {config}: FAILED ({data['error']})")
        
        # Show stabilizer information if available
        if franz_results.get('stabilizer_info'):
            print(f"  Stabilizer generator info:")
            for info in franz_results['stabilizer_info']:
                if 'error' not in info:
                    print(f"    Graph {info['graph_index'] + 1}: X={info['logical_X']}, Cost={info['cost']}")
                else:
                    print(f"    Error: {info['error']}")
                    
    else:
        print(f"  Status: FAILED")
        print(f"  Error: {franz_results['error']}")
    
    # Cost comparison
    if lx_results['success'] and franz_results['success']:
        lx_cost = lx_results['best_cost']
        franz_cost = franz_results['best_cost']
        
        print(f"\nCost Comparison:")
        print(f"  LXMixer: {lx_cost}")
        print(f"  Franz:   {franz_cost}")
        
        if lx_cost < franz_cost:
            print(f"  Winner: LXMixer (by {franz_cost - lx_cost} CNOTs)")
        elif franz_cost < lx_cost:
            print(f"  Winner: Franz (by {lx_cost - franz_cost} CNOTs)")
        else:
            print(f"  Result: TIE")
        
        # Time comparison
        print(f"\nTime Comparison:")
        print(f"  LXMixer: {lx_results['execution_time']:.4f}s")
        print(f"  Franz:   {franz_results['execution_time']:.4f}s")
        
        if lx_results['execution_time'] < franz_results['execution_time']:
            print(f"  Faster: LXMixer")
        else:
            print(f"  Faster: Franz")
        
        # Architecture comparison
        print(f"\nArchitecture Comparison:")
        print(f"  LXMixer approach: Orbit-based with stabilizer projectors")
        print(f"  Franz approach: Graph-based with commuting Pauli strings")
        print(f"  LXMixer orbits: {len(lx_results['orbits'])}")
        if franz_results.get('detailed_analyzer') and hasattr(franz_results['detailed_analyzer'], 'solution'):
            print(f"  Franz graphs: {len(franz_results['detailed_analyzer'].solution)}")
    
    print("="*80)
    
    return lx_results, franz_results

def analyze_random_samples(n, num_basis_states, num_samples=5):
    """
    Analyze multiple random samples of a given size.
    
    Args:
        n: Number of qubits
        num_basis_states: Number of basis states to include
        num_samples: Number of random samples to analyze
    """
    print(f"ANALYZING {num_samples} RANDOM SAMPLES")
    print(f"n={n} qubits, |B|={num_basis_states}")
    print("="*60)
    
    lx_costs = []
    franz_costs = []
    lx_times = []
    franz_times = []
    
    for sample_idx in range(num_samples):
        print(f"\nSAMPLE {sample_idx + 1}:")
        print("-" * 30)
        
        # Generate random sample
        B_strings, B_integers = generate_random_B(n, num_basis_states)
        print(f"B = {B_strings}")
        
        try:
            # Analyze this sample with both mixers
            lx_results, franz_results = compare_mixers(n, num_basis_states, B_strings, B_integers)
            
            if lx_results['success']:
                lx_costs.append(lx_results['best_cost'])
                lx_times.append(lx_results['execution_time'])
                print(f"LXMixer cost: {lx_results['best_cost']}")
            
            if franz_results['success']:
                franz_costs.append(franz_results['best_cost'])
                franz_times.append(franz_results['execution_time'])
                print(f"Franz cost: {franz_results['best_cost']}")
            
            if lx_results['success'] and franz_results['success']:
                improvement = franz_results['best_cost'] - lx_results['best_cost']
                print(f"Cost difference (Franz - LXMixer): {improvement}")
            
        except Exception as e:
            print(f"Error analyzing sample: {e}")
            continue
    
    if lx_costs and franz_costs:
        print(f"\nSUMMARY STATISTICS:")
        print("-" * 30)
        print(f"LXMixer costs: min={min(lx_costs)}, max={max(lx_costs)}, mean={np.mean(lx_costs):.1f}")
        print(f"Franz costs: min={min(franz_costs)}, max={max(franz_costs)}, mean={np.mean(franz_costs):.1f}")
        print(f"Average cost difference (Franz - LXMixer): {np.mean(franz_costs) - np.mean(lx_costs):.1f}")
        print(f"LXMixer times: min={min(lx_times):.3f}s, max={max(lx_times):.3f}s, mean={np.mean(lx_times):.3f}s")
        print(f"Franz times: min={min(franz_times):.3f}s, max={max(franz_times):.3f}s, mean={np.mean(franz_times):.3f}s")
        
        # Win statistics
        lx_wins = sum(1 for lx, franz in zip(lx_costs, franz_costs) if lx < franz)
        franz_wins = sum(1 for lx, franz in zip(lx_costs, franz_costs) if franz < lx)
        ties = sum(1 for lx, franz in zip(lx_costs, franz_costs) if lx == franz)
        
        print(f"Win statistics: LXMixer={lx_wins}, Franz={franz_wins}, Ties={ties}")

def interactive_mode():
    """
    Interactive mode for exploring specific cases.
    """
    print("Welcome to Single Case Mixer Comparison Tool!")
    print("=" * 50)
    
    while True:
        print("\nOptions:")
        print("1. Random B set comparison")
        print("2. Custom B set comparison")
        print("3. Predefined test cases")
        print("4. Random samples analysis")
        print("5. Exit")
        
        choice = input("\nEnter your choice (1-5): ").strip()
        
        if choice == '1':
            try:
                n = int(input("Enter number of qubits (n): "))
                m = int(input("Enter size of B set (m): "))
                
                if m > 2**n:
                    print(f"Error: m={m} cannot be larger than 2^n={2**n}")
                    continue
                    
                compare_mixers(n, m)
                
            except ValueError:
                print("Invalid input. Please enter integers.")
                
        elif choice == '2':
            try:
                n = int(input("Enter number of qubits (n): "))
                print("Enter B set as binary strings (space-separated):")
                print("Example: 000 001 010 011")
                
                b_input = input("B set: ").strip().split()
                
                # Validate input
                if not all(len(b) == n and all(c in '01' for c in b) for b in b_input):
                    print(f"Error: All strings must be {n} bits long and contain only 0s and 1s")
                    continue
                
                B_strings = b_input
                B_integers = [int(b, 2) for b in B_strings]
                m = len(B_strings)
                
                compare_mixers(n, m, B_strings, B_integers)
                
            except ValueError:
                print("Invalid input.")
                
        elif choice == '3':
            print("\nPredefined test cases:")
            print("1. n=3, B=[6,3,1,5,0,4,2] (7 elements)")
            print("2. n=4, B=[14,12,9,4,3] (5 elements, from article)")
            print("3. n=4, B=[0,15,1,13,14,12,2,3] (8 elements)")
            
            test_choice = input("Select test case (1-3): ").strip()
            
            if test_choice == '1':
                n = 3
                B_integers = [6, 3, 1, 5, 0, 4, 2]
                B_strings = [f'{b:03b}' for b in B_integers]
                m = len(B_integers)
                compare_mixers(n, m, B_strings, B_integers)
                
            elif test_choice == '2':
                n = 4
                B_integers = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
                B_strings = [f'{b:04b}' for b in B_integers]
                m = len(B_integers)
                compare_mixers(n, m, B_strings, B_integers)
                
            elif test_choice == '3':
                n = 4
                B_integers = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100, 0b0010, 0b0011]
                B_strings = [f'{b:04b}' for b in B_integers]
                m = len(B_integers)
                compare_mixers(n, m, B_strings, B_integers)
                
            else:
                print("Invalid choice.")
                
        elif choice == '4':
            try:
                n = int(input("Enter number of qubits (n): "))
                m = int(input("Enter size of B set (m): "))
                num_samples = int(input("Enter number of samples (default 5): ") or "5")
                
                if m > 2**n:
                    print(f"Error: m={m} cannot be larger than 2^n={2**n}")
                    continue
                    
                analyze_random_samples(n, m, num_samples)
                
            except ValueError:
                print("Invalid input. Please enter integers.")
                
        elif choice == '5':
            print("Goodbye!")
            break
            
        else:
            print("Invalid choice. Please enter 1-5.")

if __name__ == "__main__":
    # You can either run in interactive mode or directly call compare_mixers
    
    # Example: Direct comparison
    # n = 3
    # m = 7
    # compare_mixers(n, m)
    
    # Or predefined B set
    # n = 4
    # B_integers = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
    # B_strings = [f'{b:04b}' for b in B_integers]
    # m = len(B_integers)
    # compare_mixers(n, m, B_strings, B_integers)
    
    # Interactive mode
    interactive_mode()
