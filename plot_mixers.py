import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import numpy as np
from utils import pauli_int_to_str
import math

X0 = 0
FONTSIZE = 7

def draw_nodes(ax, x, nodes, nL):
    x_text = -0.12 + x
    N = len(nodes)
    positions = {}
    for i in range(N):
        y = N - i
        positions[i] = (0, y)
        ax.plot(x, y, 'o', markersize=16, color='white')
        ax.text(x_text, y, fr'$|${nodes[i]:0{nL}b}$\rangle$', 
                fontsize=FONTSIZE, verticalalignment='center', color='black')
        ax.set_ylim(0, N + 1)
        ax.set_xlim(-1, 1)
        ax.axis('off')
        
def draw_arc(ax, start, end, r=0.2, color='black', lw=1):
    """
    Draw an arc from start to end with a given radius.
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
        
def draw_orbit(ax, nodes, Xs, LX, x, color, r, lw):
    for X in Xs:
        for node in nodes:
            for neighbor in nodes:
                if node < neighbor and (node, neighbor) in LX.family_of_valid_graphs[X]:
                    if node%2 == 0:
                        y0 = LX.nB - node
                        y1 = LX.nB - neighbor
                    else:
                        y0 = LX.nB - neighbor
                        y1 = LX.nB - node
                    start = (x, y0)
                    end = (x, y1)
                    draw_arc(ax, start, end, r=r, color=color, lw=lw)

def get_colors(n, cmap_name="tab20", start=0.0, end=1.0):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i) for i in np.linspace(start, end, n)]
                    
def draw_mixer_graph(ax, combination, Xs, LX, x=X0, r=0.3, lw=1, cmap="tab20"):
    draw_nodes(ax, x, LX.B, LX.nL)

    colors = get_colors(len(combination), cmap_name=cmap)
        
    orbit_labels = []
    for color_n, nodes in enumerate(combination):
        orbit_labels.append(fr"$\langle${", ".join(f"${pauli_int_to_str(X, LX.nL)}$"for X in Xs[color_n])}$\rangle$")
        draw_orbit(ax, nodes, Xs[color_n], LX, x, colors[color_n], r+color_n*0.003, lw)
                    
    handles = [plt.Line2D([0], [0],
                    linestyle="-",      # solid line
                    color=color,        # line color
                    lw=1,               # line width
                    label=label)
        for color, label in zip(colors, orbit_labels)]
    ax.legend(
        handles=handles,
        fontsize=FONTSIZE,
        title_fontsize=FONTSIZE,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),  # Center above the plot in axes coords
        bbox_transform=ax.transAxes,  # Use axes coords, not data coords
        frameon=False
    )
    
def draw_best_graphs(LX, x=X0, r=0.3, lw=1, cmap="tab20"):
    N_plots = len(LX.best_combinations)
    # fig, ax = plt.subplots(1, N_plots, figsize=(int(math.log2(LX.nB))*N_plots, LX.nB+max([len(graph_Xs) for Xs in LX.best_Xs for graph_Xs in Xs])*0.5))
    fig, ax = plt.subplots(1, N_plots, figsize=(int(math.log2(LX.nB))*N_plots, LX.nB))

    if not isinstance(ax, (list, np.ndarray)):  # Ensure ax is iterable
        ax = [ax]
    
    for plot_n, combination in enumerate(LX.best_combinations):
        draw_mixer_graph(ax[plot_n], combination, LX.best_Xs[plot_n], LX, x=x, r=r, lw=lw, cmap=cmap)
        # ax[plot_n].set_ylim(0, LX.nB + legend_height)

def group_family_of_valid_graphs(family_of_valid_graphs, group_size):
    """
    Groups the family of valid graphs into smaller groups based on the specified group size.
    """
    grouped_graphs = []
    current_group = {}
    for i, (X, edges) in enumerate(family_of_valid_graphs.items()):
        current_group[X] = edges
        if (i + 1) % group_size == 0 or i == len(family_of_valid_graphs) - 1:
            grouped_graphs.append(current_group)
            current_group = {}
    return grouped_graphs

def draw_family_of_valid_graphs(LX, x=X0, r=0.3, lw=1, group_size=3, cmap="tab20"):
    """
    Draws the family of valid graphs for the logical X operators, grouped into fewer graphs.
    """
    grouped_graphs = group_family_of_valid_graphs(LX.family_of_valid_graphs, group_size)
    N_plots = len(grouped_graphs)
    # fig, ax = plt.subplots(1, N_plots, figsize=(N_plots * 3, LX.nB+group_size*0.5))
    fig, ax = plt.subplots(1, N_plots, figsize=(int(math.log2(LX.nB))*N_plots, LX.nB))
    
    if not isinstance(ax, (list, np.ndarray)):  # Ensure ax is iterable
        ax = [ax]
    colors = get_colors(len(LX.family_of_valid_graphs), cmap_name=cmap, start=0.0, end=1.0)
    X_color = 0
    for plot_n, group in enumerate(grouped_graphs):
        group_colors = []
        for X, edges in group.items():
            group_colors.append(colors[X_color])
            draw_nodes(ax[plot_n], X0, LX.B, LX.nL)
            for edge in edges:
                if edge[0] % 2 == 0:
                    y0 = LX.nB - edge[0]
                    y1 = LX.nB - edge[1]
                else:
                    y0 = LX.nB - edge[1]
                    y1 = LX.nB - edge[0]
                start = (X0, y0)
                end = (X0, y1)
                draw_arc(ax[plot_n], start, end, color=colors[X_color], r=r, lw=lw)
            X_color += 1
        
        group_labels = [pauli_int_to_str(X, LX.nL) for X in group.keys()]
        handles = [plt.Line2D([0], [0],
                linestyle="-",      # solid line
                color=color,        # line color
                lw=1,               # line width
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