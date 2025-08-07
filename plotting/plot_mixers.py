import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import numpy as np
from utils import pauli_int_to_str
import math

X0 = 0
FONTSIZE = 7

class Plotter:
    """
    A class to handle drawing of logical X mixers and their orbits.
    """
    def __init__(self, LX):
        """
        Initializes the Plotter with a logical X mixer instance.

        Args:
            LX (LXMixer): An instance of the LXMixer class containing `B`, `nB`, `nL`, `family_of_valid_graphs`, `orbits` (with projectors and costs)
            and best solutions.
        """
        self.LX = LX
        self.x0 = X0
        self.fontsize = FONTSIZE

    def draw_nodes(self, ax):
        """
        Draws nodes on the plot at a given x position.

        Args:
            ax (matplotlib.axes.Axes): The axes on which to draw the nodes.
        """
        for i in range(self.LX.nB):
            y = self.LX.nB - i # y-coordinate is inverted to match the plot's y-axis direction.
            ax.plot(self.x0, y, "o", markersize=16, color="white") # Draw the node as a white circle.
            ax.text(self.x0, y, fr"$|${self.LX.B[i]:0{self.LX.nL}b}$\rangle$", 
                    fontsize=FONTSIZE, verticalalignment="center", horizontalalignment="center", color="black")
            ax.set_xlim(-1, 1)
            ax.axis("off")  # Hide the axes for a cleaner look.
        
    def draw_arc(self, ax, start, end, r=0.2, color="black", lw=1.0):
        """
        Draw an arc from start to end with a given radius.
        Args:
            ax (matplotlib.axes.Axes): The axes on which to draw the arc.
            start (tuple): Starting point of the arc (x, y).
            end (tuple): Ending point of the arc (x, y).
            r (float): Radius of the arc.
            color (str): Color of the arc.
            lw (float): Line width of the arc.
        """
        x0, y0 = start
        x1, y1 = end
        dx = x1 - x0
        dy = y1 - y0
        
        ctrl_x = (x0 + x1)/2 - r*dy
        ctrl_y = (y0 + y1)/2 + r*dx
        
        path = Path([start, (ctrl_x, ctrl_y), end], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        patch = FancyArrowPatch(
            path=path,
            arrowstyle="-",  # Default is undirected, no arrowheads
            color=color,
            lw=lw,
            connectionstyle="arc3",
            # alpha=0.8
        )
        ax.add_patch(patch)
        
    def draw_orbit(self, ax, nodes, Xs, color, r, lw):
        """
        Draws the orbit of nodes connected by logical X operators.

        Args:
            ax (matplotlib.axes.Axes): The axes on which to draw the orbit.
            nodes (set): Set of nodes in the orbit.
            Xs (list): List of logical X operators (int representations) that connect the nodes.
            color (str): Color for the arcs connecting the nodes.
            r (float): Radius for the arcs.
            lw (float): Line width for the arcs.
        """
        for X in Xs:
            for node in nodes:
                side = 0
                for neighbor in nodes:
                    side += 1
                    if node < neighbor and (node, neighbor) in self.LX.family_of_valid_graphs[X]:
                        if side%2 == 0:
                            y0 = self.LX.nB - node
                            y1 = self.LX.nB - neighbor
                        else:
                            y0 = self.LX.nB - neighbor
                            y1 = self.LX.nB - node
                        start = (self.x0, y0)
                        end = (self.x0, y1)
                        self.draw_arc(ax, start, end, r=r, color=color, lw=lw)

    def get_colors(self, n, cmap_name="tab20", start=0.0, end=1.0):
        """
        Generate a list of colors from a colormap.
        
        Args:
            n (int): Number of colors to generate.
            cmap_name (str): Name of the colormap to use.
            start (float): Start point in the colormap (0.0 to 1.0).
            end (float): End point in the colormap (0.0 to 1.0).
        
        Returns:
            List of colors.
        """
        cmap = plt.get_cmap(cmap_name)
        return [cmap(i) for i in np.linspace(start, end, n)]
                    
    def draw_mixer_graph(self, ax, combination, Xs, r=0.3, lw=1, cmap="tab20"):
        """
        Draws the mixer graph for a given combination of nodes.
        
        Args:
            ax (matplotlib.axes.Axes): The axes on which to draw the graph.
            combination (list): List of sets, where each set contains nodes in an orbit.
            Xs (list): List of logical X operators (int representations) that connect the nodes.
            r (float): Radius for the arcs.
            lw (float): Line width for the arcs.
            cmap (str): Name of the colormap to use for coloring the orbits.
        """
        
        self.draw_nodes(ax) # Draw the nodes of B.

        colors = self.get_colors(len(combination), cmap_name=cmap) # Get colors for each orbit.
            
        orbit_labels = [] # Labels for the legend.
        for color_n, nodes in enumerate(combination):
            orbit_labels.append(fr"$\langle${", ".join(f"${pauli_int_to_str(X, self.LX.nL)}$"for X in Xs[color_n])}$\rangle$") # Orbit label <X_1, ..., X_l>.
            self.draw_orbit(ax, nodes, Xs[color_n], colors[color_n], r+color_n*0.005, lw) # Draw the orbit with the corresponding color.
                        
        handles = [plt.Line2D([0], [0],
                        linestyle="-",
                        color=color,
                        lw=lw,
                        label=label)
            for color, label in zip(colors, orbit_labels)]
        ax.legend(
            handles=handles,
            fontsize=FONTSIZE,
            title_fontsize=FONTSIZE,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.02),  # Center above the plot in axes coords.
            bbox_transform=ax.transAxes,  # Use axes coords, not data coords.
            frameon=False
        )
    
    def draw_best_graphs(self, x=X0, r=0.3, lw=1, cmap="tab20", saveas=None):
        N_plots = len(self.LX.best_combinations)
        fig, ax = plt.subplots(1, N_plots, figsize=(int(math.log2(self.LX.nB))*N_plots, self.LX.nB))
        if not isinstance(ax, (list, np.ndarray)):  # Ensure ax is iterable, even if there's only one subplot/best solution.
            ax = [ax]
        for plot_n, combination in enumerate(self.LX.best_combinations):
            self.draw_mixer_graph(ax[plot_n], combination, self.LX.best_Xs[plot_n], r=r, lw=lw, cmap=cmap) # Draw the mixer graph for each best combination.

        plt.tight_layout()
        if saveas: plt.savefig(saveas)

    def group_family_of_valid_graphs(self, group_size):
        """
        Groups the family of valid graphs into smaller groups based on the specified group size.
        
        Args:
            family_of_valid_graphs (Dict[int, List[Tuple[int,...]]]): A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator.
            group_size (int): The number of logical X operators to include in each group.
            
        Returns:
            List[Dict[int, List[Tuple[int,...]]]]: A list of dictionaries, where each dictionary contains a group of logical X operators and their corresponding edges.
        """
        grouped_graphs = []
        current_group = {}
        for i, (X, edges) in enumerate(self.LX.family_of_valid_graphs.items()):
            current_group[X] = edges
            if (i + 1) % group_size == 0 or i == len(self.LX.family_of_valid_graphs) - 1:
                grouped_graphs.append(current_group)
                current_group = {}
        return grouped_graphs

    def draw_family_of_valid_graphs(self, r=0.3, lw=1, group_size=3, cmap="tab20", saveas=None):
        """
        Draws the family of valid graphs for the logical X operators, grouped into fewer graphs.
        
        Args:
            LX: An instance of the logical X mixer class containing `family_of_valid_graphs`
            x (float): The x-coordinate for the nodes.
            r (float): Radius for the arcs.
            lw (float): Line width for the arcs.
            group_size (int): The number of logical X operators to include in each group.
            cmap (str): Name of the colormap to use for coloring the orbits.
        """
        grouped_graphs = self.group_family_of_valid_graphs(group_size)
        N_plots = len(grouped_graphs)
        fig, ax = plt.subplots(1, N_plots, figsize=(int(math.log2(self.LX.nB))*N_plots, self.LX.nB))
        
        if not isinstance(ax, (list, np.ndarray)):  # Ensure ax is iterable
            ax = [ax]
        colors = self.get_colors(len(self.LX.family_of_valid_graphs), cmap_name=cmap, start=0.0, end=1.0)
        X_color = 0
        for plot_n, group in enumerate(grouped_graphs):
            group_colors = []
            for X, edges in group.items():
                group_colors.append(colors[X_color])
                self.draw_nodes(ax[plot_n])
                for edge in edges:
                    if edge[0] % 2 == 0:
                        y0 = self.LX.nB - edge[0]
                        y1 = self.LX.nB - edge[1]
                    else:
                        y0 = self.LX.nB - edge[1]
                        y1 = self.LX.nB - edge[0]
                    start = (X0, y0)
                    end = (X0, y1)
                    self.draw_arc(ax[plot_n], start, end, color=colors[X_color], r=r, lw=lw)
                X_color += 1
            
            group_labels = [pauli_int_to_str(X, self.LX.nL) for X in group.keys()]
            handles = [plt.Line2D([0], [0],
                    linestyle="-",
                    color=color,
                    lw=1,
                    label=label)
                for color, label in zip(group_colors, group_labels)]
            ax[plot_n].legend(
                handles=handles,
                fontsize=FONTSIZE,
                title_fontsize=FONTSIZE,
                loc="lower center",
                bbox_to_anchor=(0.5, 1.02),  # Center above the plot in axes coords
                bbox_transform=ax[plot_n].transAxes,  # Use the specific Axes object
                frameon=False
            )
            
        plt.tight_layout()
        if saveas: plt.savefig(saveas)