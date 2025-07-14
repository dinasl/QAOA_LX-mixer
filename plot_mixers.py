import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import numpy as np

# from Mixer import LXMixer

def plot_mixer_graph(LX):
    
    x = -0.1
    x_text = -0.15 + x
    N_plots = len(LX.best_combinations)
    fig, ax = plt.subplots(1, N_plots, figsize=(2*N_plots, LX.nB))
    
    positions = {}
    
    for plot_n, combination in enumerate(LX.best_combinations):
        for i in range(LX.nB):
            
            # Plot all nodes in B
            y = LX.nB - i
            positions[i] = (0, y)
            ax[plot_n].plot(x, y, 'o', markersize=15, color='white')
            ax[plot_n].text(x_text, y, fr'$|${LX.B[i]:0{LX.nL}b}$\rangle$', fontsize=6, verticalalignment='center', color='black')
            ax[plot_n].set_ylim(0, LX.nB + 1)
            ax[plot_n].set_xlim(-1, 1)
            ax[plot_n].axis('off')

        colors = get_colors(len(combination), cmap_name="tab20")
        for color_n, nodes in enumerate(combination):
            orbit_label = fr"$\langle${", ".join(f"{X:0{LX.nL}b}"for X in LX.orbits[nodes].Xs)}$\rangle$"
            for node in nodes:
                y0 = LX.nB - node
                for neighbor in nodes:
                    if neighbor != node:
                        y1 = LX.nB - neighbor
                        start = (x, y0)
                        end = (x, y1)
                        draw_arc(ax[plot_n], start, end, orbit_label, r=0.3, color= colors[color_n])
    
    plt.tight_layout()
    plt.show()

def draw_arc(ax, start, end, orbit_label, r=0.2, color='black', lw=1):
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
    # ax.text(ctrl_x, ctrl_y, orbit_label, fontsize=6, verticalalignment='center', horizontalalignment='center', color=color)
    
def get_colors(n, cmap_name="tab20", start=0.0, end=1.0):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i) for i in np.linspace(start, end, n)]