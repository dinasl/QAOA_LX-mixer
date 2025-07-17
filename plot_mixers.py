import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import numpy as np
import math
from utils import pauli_int_to_str

def draw_nodes(ax, x, nodes, nL):
    x_text = -0.15 + x
    N = len(nodes)
    positions = {}
    for i in range(N):
        y = N - i
        positions[i] = (0, y)
        ax.plot(x, y, 'o', markersize=15, color='white')
        ax.text(x_text, y, fr'$|${nodes[i]:0{nL}b}$\rangle$', 
                fontsize=6, verticalalignment='center', color='black')
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
                    
def draw_mixer_graph(ax, combination, Xs, LX, x, r=0.3, lw=1):
    draw_nodes(ax, x, LX.B, LX.nL)

    colors = get_colors(len(combination), cmap_name="tab20")
    # radii = get_radii(len(combination), r_min=0.1, r_max=0.3)
    
    orbit_labels = []
    for color_n, nodes in enumerate(combination):
        orbit_labels.append(fr"$\langle${", ".join(f"${pauli_int_to_str(X, LX.nL)}$"for X in Xs[color_n])}$\rangle$")
        # if color_n % 2 == 0:
        #     side = "left"
        # else:
        #     side = "right"
        draw_orbit(ax, nodes, Xs[color_n], LX, x, colors[color_n], r, lw)
                    
    handles = [plt.Line2D([0], [0],
                    linestyle="-",      # solid line
                    color=color,        # line color
                    lw=1,               # line width
                    label=label)
        for color, label in zip(colors, orbit_labels)]
    ax.legend(handles=handles, title = "Orbits", fontsize=6, title_fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.05), frameon=False)

def draw_best_graphs(LX, r=0.3, lw=1):
    x = -0.1
    N_plots = len(LX.best_combinations)
    fig, ax = plt.subplots(1, N_plots, figsize=(2*N_plots, LX.nB))
    if not isinstance(ax, (list, np.ndarray)):  # Ensure ax is iterable
        ax = [ax]
    
    for plot_n, combination in enumerate(LX.best_combinations):
        draw_mixer_graph(ax[plot_n], combination, LX.best_Xs[plot_n], LX, x, r, lw)

    plt.tight_layout()
    plt.show()
    
# def get_radii(n, r_min=0.1, r_max=0.3):
#     """
#     Generate a list of radii for arcs.
#     """
#     return [r_min + (r_max - r_min) * i / (n - 1) for i in range(n)]