# Visualization

import os, glob
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch
import seaborn as sns

import rasterio
import shapely
import geopandas as gpd
from rasterio import features
from shapely.geometry import Polygon
from skimage import measure

# Clip image rasters for enhanced visualization
def clip(img, percentile):
    out = np.zeros_like(img.shape[2])
    for i in range(img.shape[2]):
        a = 0 
        b = 255 
        c = np.percentile(img[:,:,i], percentile)
        d = np.percentile(img[:,:,i], 100 - percentile)        
        t = a + (img[:,:,i] - c) * (b - a) / (d - c)    
        t[t<a] = a
        t[t>b] = b
        img[:,:,i] =t
    rgb_out = img.astype(np.uint8)   
    return rgb_out

# Read ASCII
def read_asc_file(file_path):
    with open(file_path, 'r') as file:
        # Skip the header lines (typically six lines in an ASC file)
        for _ in range(6):
            next(file)
        
        # Read the data as a 2D list
        data = []
        for line in file:
            row_data = list(map(float, line.strip().split()))
            data.append(row_data)
    
    return data

def asc_to_numpy(asc_data):
    return np.array(asc_data)

# Single FARSITE path directories
def load_farsite(farsite_path):

    # Open Farsite path
    farsite_perim = gpd.read_file(farsite_path)
    farsite_perim.crs = None

    return farsite_perim

def farsite_contour(farsite_img):

    # Find the contours in the binary mask
    contours = measure.find_contours(farsite_img, 0.8)

    # Increase the number of points on the contour
    increased_contours = []
    for contour in contours:
        increased_contour = measure.approximate_polygon(contour, tolerance=0)
        increased_contours.append(increased_contour)

    return increased_contours

# Multiple FARSITE path directories (multiple time intervals)
def load_farsite_list(farsite_path_list):

    data, data_stack = [], []

    # Loop through all paths (based on number of hours burned)
    for farsite_path in farsite_path_list:
        data = pd.read_csv(farsite_path, header=None, delimiter=' ')
        # Convert the data to a NumPy array
        data = np.array(data)
        data_stack.append(data) # List stacked with arrays

    return data_stack

# Load Cell2Fire csv results
def load_cell2fire(cell2fire_path):

    # Open Cell2Fire CSV data
    data = []
    with open(cell2fire_path, 'r') as file:
        for line in file:
            values = line.strip().split(' ')
            data.append([int(value) for value in values])

    # Convert the data to a NumPy array
    data = np.array(data)

    return data
    
def rasterize_cell2fire(polygon, size):
            
    burn_value = 0
    xmin, ymin, xmax, ymax = (0, 0, size, size)
    binary_grid = np.zeros((size, size), dtype=np.uint8)

    rasterized = features.rasterize([(polygon.geometry, burn_value)], out_shape=(size, size), fill=1) # Inverted
    binary_grid = np.where(rasterized == 1, 1, binary_grid)

    return binary_grid



def compare_single_grids_v2(base_path, grid1, grid2, size, instance, save=False):

    """
    Compare fire spread OVER TIME between two grids and visualize the difference.

    Parameters:
    - base_path (str): Base directory path for saving images.
    - grid1 (list): List of Cell2FireML's 2D grids
    - grid2 (list): List of FarSite's 2D grids
    - size (int): Size of the grid.
    - instance (str): Name of instance (e.g., "f101_100_ws10_CS100_5HR" which means:
                                        Fuel type 101, grid size 100, wind speed 10, 
                                        cell size = 100, fire duration = 5 hours)
    - save (bool): Flag indicating whether to save the plot as an image (Default=False).
    """

    # Set colors for fire progression over time
    colors = ['#579c8c', '#b5e6b5', '#f7f7f7', '#ffc872', '#ff9696']    
    cmap1 = ListedColormap([colors[num] for num in range(len(colors))]) # Colors for plot 1 (fire spread progression)
    cmap2 = ListedColormap(['firebrick','white','royalblue'])           # Colors for plot 2 (over/under estimation)

    # Plot
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 12), sharex=True, sharey=True)

    ### Plot 1: Comparison of fire spread ###
    # Plot Cell2Fire
    axes[0].set_xlim([0,size]); axes[0].set_ylim([0,size])

    axes[0].imshow(grid1[4], cmap=plt.cm.colors.ListedColormap(['none', colors[4]]), interpolation='nearest', vmin=0, vmax=1, label= '1HR')
    axes[0].imshow(grid1[3], cmap=plt.cm.colors.ListedColormap(['none', colors[3]]), interpolation='nearest', vmin=0, vmax=1, label= '2HR')
    axes[0].imshow(grid1[2], cmap=plt.cm.colors.ListedColormap(['none', colors[2]]), interpolation='nearest', vmin=0, vmax=1, label= '3HR')
    axes[0].imshow(grid1[1], cmap=plt.cm.colors.ListedColormap(['none', colors[1]]), interpolation='nearest', vmin=0, vmax=1, label= '4HR')
    axes[0].imshow(grid1[0], cmap=plt.cm.colors.ListedColormap(['none', colors[0]]), interpolation='nearest', vmin=0, vmax=1, label= '5HR')

    # Plot FARSITE
    for num, grid in enumerate(grid2):
        contours = farsite_contour(grid)
        for i, contour in enumerate(contours):
            axes[0].plot(contour[:, 1], contour[:, 0], linewidth=1, color='k')

    # Set extents, scale bar, and colorbar
    scalebar = ScaleBar(1, 'm', length_fraction=0.2, width_fraction=0.02, frameon=False, location='lower right', label='Scale', font_properties={"size": 20})
    axes[0].add_artist(scalebar)

    divider1 = make_axes_locatable(axes[0])
    cax1 = divider1.append_axes('bottom', size='5%', pad=0.1)
    sm1 = plt.cm.ScalarMappable(cmap=cmap1)
    cbar1 = fig.colorbar(sm1, cax=cax1, orientation='horizontal')
    cbar1.set_ticks([0.1, 0.3, 0.5, 0.7, 0.9])
    cbar1.set_ticklabels(['1H','2H','3H','4H','5H'], size=18)

    ### Plot 2: Difference in fire spread ###
    diff = grid2[-1] - grid1[-1]
    axes[1].imshow(diff, cmap=cmap2, interpolation='nearest', vmin=-1, vmax=1)

    # Set extents, scale bar, and colorbar
    axes[1].set_xlim([0,size]); axes[1].set_ylim([0,size])
    scalebar = ScaleBar(1, 'm', length_fraction=0.2, width_fraction=0.02, frameon=False, location='lower right', label='Scale', font_properties={"size": 20})
    axes[1].add_artist(scalebar)

    divider2 = make_axes_locatable(axes[1])
    cax2 = divider2.append_axes('bottom', size='5%', pad=0.1)
    sm2 = plt.cm.ScalarMappable(cmap=cmap2)
    sm2.set_clim(-1, 1)  # Set the color limits to match the desired values
    cbar2 = fig.colorbar(sm2, cax=cax2, orientation='horizontal')
    cbar2.set_ticks([2/3, 0, -2/3])
    cbar2.set_ticklabels(['Underestimate','No Difference','Overestimate'],size=16, ha='center')  # Set custom tick labels
    cbar2.ax.yaxis.set_tick_params(pad=20)

    # Adjust grid extents 
    if size == 500:
        axes[0].set_xlim([0,200]); axes[0].set_ylim([150,350])
        axes[1].set_xlim([0,200]); axes[1].set_ylim([150,350])

    axes[0].set_xticks([]);axes[0].set_yticks([])
    axes[1].set_xticks([]);axes[1].set_yticks([])

    fig.tight_layout()

    if save:
        plt.savefig(os.path.join(base_path, instance) + '_solid_v2.png', bbox_inches='tight', dpi=300)
    plt.show()

        

def compare_single_grids(base_path, grid1, grid2, x_size, y_size, instance, fuel_type=None, wind_field=False, save=False):

    """
    Compare fire spread between two grids and visualize the difference for ONLY the final burned time.

    Parameters:
    - base_path (str): Base directory path for saving images.
    - grid1 (np.array): Cell2FireML's 2D grid
    - grid2 (list): List of FarSite's 2D grids
    - size (int): Size of the grid.
    - instance (str): Name of instance (e.g., "f101_100_ws10_CS100_5HR" which means:
                                        Fuel type 101, grid size 100, wind speed 10, 
                                        cell size = 100, fire duration = 5 hours)
    - save (bool): Flag indicating whether to save the plot as an image (Default=False).
    """

    # Set colors for fire progression
    colors = ['#579c8c', '#b5e6b5', '#f7f7f7', '#ffc872', '#ff9696']    # Colors for plot 1 (fire spread progression)
    cmap2 = ListedColormap(['firebrick','white','royalblue'])           # Colors for plot 2 (over/under estimation)

    # Plot
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 15))
    
    # Adjust subplots to make them equal sizes
    fig.subplots_adjust(wspace=0.42)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.05])
    axes = [plt.subplot(gs[0]), plt.subplot(gs[1])]

    if wind_field:
    
        # Plot the wind field using quiver plot
        x = np.linspace(0, x_size, int(x_size/16))
        y = np.linspace(0, y_size, int(y_size/16))
        X, Y = np.meshgrid(x, y)
        U = np.ones_like(X)  # Constant eastward wind component
        V = np.zeros_like(Y)  # No northward wind component
        axes[0].quiver(X, Y, U, V, color='grey', alpha=0.3, zorder=0)

        # Wind direction subplot
        inner_ax = axes[0].inset_axes([0.85, 0.9, 0.15, 0.1])
        inner_ax.annotate("",
                    xy=(0.85,0.35),  # starting point (in axis coordinates)
                    xytext=(0.25, 0.35),  # ending point (in axis coordinates)
                    arrowprops=dict(arrowstyle="->", color="black", linewidth=2))
        inner_ax.annotate("Wind", xy=(0.46, 0.55), fontsize=18, ha='center')
        inner_ax.set_xticklabels([]); inner_ax.set_yticklabels([])

    ### Plot 1: Comparison of fire spread ###
    # Plot Cell2Fire
    axes[0].set_xlim([0,x_size]); axes[0].set_ylim([y_size,0])
    axes[0].imshow(grid1, cmap=plt.cm.colors.ListedColormap(['none', colors[-1]]), interpolation='nearest', vmin=0, vmax=1, label= '5HR')

    # Plot FARSITE
    for num, grid in enumerate(grid2):
        contours = farsite_contour(grid) # Compute FARSITE perimeters
        for i, contour in enumerate(contours):
            axes[0].plot(contour[:, 1], contour[:, 0], linewidth=1, color='k')

    # Scalebar
    scalebar = ScaleBar(1, 'm', length_fraction=0.2, height_fraction=0.02, frameon=False, location='lower right', label='Scale', font_properties={"size": 18})
    axes[0].add_artist(scalebar)

    ### Plot 2: Difference in fire spread ###
    diff = grid2[-1] - grid1
    axes[1].imshow(diff, cmap=cmap2, interpolation='nearest', vmin=-1, vmax=1)

    # Set extents, scale bar, and colorbar
    axes[1].set_xlim([0,x_size]); axes[1].set_ylim([y_size,0])
    scalebar = ScaleBar(1, 'm', length_fraction=0.2, width_fraction=0.02, frameon=False, location='lower right', label='Scale', font_properties={"size": 18})
    axes[1].add_artist(scalebar)

    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    sm = plt.cm.ScalarMappable(cmap=cmap2)
    sm.set_clim(-1, 1)  # Set the color limits to match the desired values
    cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
    cbar.set_ticks([2/3, 0, -2/3])
    cbar.set_ticklabels(['Underestimation','No Difference','Overestimation'], size=14, rotation=90,  va='center')  # Set custom tick labels

    if fuel_type == 101 and x_size == 500:
        axes[0].set_xlim([0,200]); axes[0].set_ylim([150,350])
        axes[1].set_xlim([0,200]); axes[1].set_ylim([150,350])
    elif fuel_type == 102:
        axes[0].set_xlim([0,500]); axes[0].set_ylim([0,500])
        axes[1].set_xlim([0,500]); axes[1].set_ylim([0,500])

    fig.tight_layout()

    if save:
        plt.savefig(os.path.join(base_path, instance) + '_solid_v2.png', bbox_inches='tight', dpi=300)
    plt.show()



