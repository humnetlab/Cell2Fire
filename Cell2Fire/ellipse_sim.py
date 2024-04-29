## Step 2.1 Simulate reference image (theoretical)
import networkx as nx
import numpy as np
import os, glob, tqdm
from Cell2Fire.ellipse_algorithms import EllipseAlgorithm

# Creates a simple 8 neighbors grid graph with nxm cells per side and cellside [m] 
def create_graph(nxm=100, cellside=100):
    # Basic G
    G = nx.grid_2d_graph(nxm, nxm)  # 4x4 grid
    
    # Add edges and set all weights to 1
    for edge in G.edges:
        G.edges[edge]['weight'] = cellside
        G.edges[edge]['len'] = cellside

    # Add diagonals
    G.add_edges_from([
        ((x, y), (x + 1, y + 1))
        for x in range(nxm - 1)
        for y in range(nxm - 1)
    ] + [
        ((x+1, y), (x, y + 1))
        for x in range(nxm - 1)
        for y in range(nxm - 1)
    ], weight=np.sqrt(2) * cellside)
    
    # Basic layour
    pos = dict( (n, n) for n in G.nodes() )
    labels = dict( ((i, j), i * nxm + j) for i, j in G.nodes() )

    # Return graph
    return G, pos, labels 

## Elliptical fittings
# Alexander elliptical definition using HROS m/min and WS in km/hr
def ellipse_definition_alexander(ros_mmin=0., ws_kmhr=0., verbose=False):
    # b: semi-major | a: semi-minor | c: ignition to center (centre to foci)
    ws_ms = ws_kmhr * 1000 / 3600
    if verbose:
        print('ws [Km/h]:', ws_kmhr)
        print('ws [m/s]:', ws_ms)
        print('ros [m/min]:', ros_mmin)

    # Fit ellipse
    a, b, c = EllipseAlgorithm.alexander(ros=ros_mmin, wind=ws_ms)

    # Info
    if verbose:
        print('Semi-major (b):', b)
        print('Semi-minor (a):', a)
        print('c:', c)
        
    # Return values
    return a, b, c

# Basic elliptical distribution using provided HROS m/min and LB
def ellipse_definition_tester(ros_mmin=0., lb=0., verbose=False):
    # b: semi-major | a: semi-minor | c: ignition to center (centre to foci)
    if verbose:
        print('lb:', lb)
        print('ros [m/min]:', ros_mmin)

    # Fit ellipse
    a, b, c = EllipseAlgorithm.tester(ros=ros_mmin, lb=lb)

    # Info
    if verbose:
        print('Semi-major (b):', b)
        print('Semi-minor (a):', a)
        print('c:', c)
        
    # Return values
    return a, b, c

# FBP ellipse using HROS, FROS and BROS in the same unit [e.g., m/min]
def ellipse_fbp(forward, flank, back, verbose=False):
    # X, Y 
    print("Initial Data to fit ellipse")
    a = (forward + back) / 2
    _x = [0.0, 
          back,
          back,
          (forward + back) / 2., 
          (forward + back) / 2., 
          (forward + back)]
    _y = [0.0, 
          np.power(flank, 2) / a,
          - (np.power(flank, 2) / a), 
          flank, 
          -flank, 
          0.0]

    # To numpy
    _x = np.array(_x)
    _y = np.array(_y)
    
    # Fit the Ellipse
    if verbose:
        print('Fit ellipse')
    FBP_Ellipse = Ellipse_FBP(_x, _y, verbose=False)

    if verbose:
        print('Get parameters')
    a, b = FBP_Ellipse.get_parameters();
    # print('Semi-major:', a)
    # print('Semi-minor:', b)

    # Placeholder for c (we don't need it)
    return a, b, None


## Elliptical distributions
# Parameterization from foci: a = semi-major | b = semi-minor 
def rhoTheta(theta, a, b):
    """Ellipse with ignition cell at one of the foci"""
    # Pi value
    pi = 3.141592653589793
    
    # Calculate c
    c2 = a ** 2 - b ** 2
    
    # Eccentricity
    e = np.sqrt(c2) / a

    # Effective Rho 
    r1 = a * (1 - e ** 2);
    r2 = 1 - e * np.cos(theta * pi / 180.0)
    r = r1 / r2;
    
    # return Rho
    return r

# Parameterization from center: a: semi-minor | b: semi-major
def rhoTheta_center(theta, a, b):
    """Ellipse centered in the ignition cell"""
    theta_rads = np.deg2rad(theta)
    r = (a * b) / np.sqrt((b * np.sin(theta_rads)) ** 2 + (a * np.cos(theta_rads)) ** 2)
    
    # Return Rho
    return r

# Parameterization from foci in radians: a: semi-major | b: semi-minor
def rhoTheta_foci(theta, a, b):
    """Ellipse centered in the foci (ignition) cell [equivalent to rhoTheta fn in radians]"""
    theta_rads = np.deg2rad(theta)
    e = np.sqrt(1 - (b ** 2)/(a ** 2))
    r = ( a * (1 - e ** 2) ) / ( 1 -  e * np.cos(theta_rads) )
    
    # Return Rho
    return r


## Angle calculations
# Returns angle in degrees from u to v (nodes in the graph according to the relabel and number of nodes)
def node2angle(u, v, factor=100):
    angle = None
    if v == (u + factor):
        angle = 0
    elif v == (u + factor + 1):
        angle = 45
    elif v == (u + 1):
        angle = 90
    elif v == (u - factor + 1):
        angle = 135
    elif v == (u - factor):
        angle = 180
    elif v == (u - factor - 1):
        angle = 225
    elif v == (u - 1):
        angle = 270
    elif v == (u + factor - 1):
        angle = 315
    
    if angle is None:
        print(u,'to', v)
    
    return angle


## Angle-to-ros distribution wrapper
# Map angle to ROS value
def angle2ros_fn(degree, thetafire, a, b, centered_at_cell=False, verbose=False):
    # Adjust ellipse for thetafire
    degree -= thetafire
    
    # Adjust
    if degree >= 360:
        degree -= 360
    if degree < 0:
        degree += 360
        
    # Foci parameterization
    if centered_at_cell is False:
        r = rhoTheta(degree, b, a) 
    # Center parameterization
    else:
        r = rhoTheta_center(degree, a, b)
    
    # Debug
    if degree in [0, 45, 90, 135, 180, 225, 270, 315, 360] and verbose:    
        print('pho(', degree, '):', r)
        
    return degree, r

# Auxiliary function for fast debugging out of C++ of different distribution schemes
def simulator_for_debugging(G, 
                            n_ignition=5050,
                            fperiods=60,
                            fperiod_length=1,
                            lb=None,
                            hros_mmin=18.94,
                            fros_mmin=None,
                            bros_mmin=None,
                            ws_kmhr=10,
                            theta_head_fire=90.,
                            ellipse_fitting='alexander',
                            verbose=True,
                            print_ros_values=False,
                            nxm=100,
                            labels=None
                           ):

    """Simulator for debugging with constant weather conditions per time step"""
    # Basic info
    # print('-' * 20, 'Fire Simulator', '-' * 20)
    # print('Ignition cell:', n_ignition)
    # print('Fperiods:', fperiods)
    # print('Fperiods length [min]:', fperiod_length)
    # print('Elliptical scheme:', ellipse_fitting)
    
    # Fire progress
    fprogress = {_:0 for _ in G.nodes()}
    onfire = []
    n_colors = {_:'green' for _ in G.nodes()}

    # Params
    thetafire = theta_head_fire
    centered_at_cell = True if 'centered' in ellipse_fitting else False

    # Ellipse fitting
    if ellipse_fitting.lower() in ['alexander', 'alexander_centered']:
        # Ellipse fitting using HROS and wind speed in km/hr
        a, b, c = ellipse_definition_alexander(ros_mmin=hros_mmin,
                                               ws_kmhr=ws_kmhr,
                                               verbose=verbose)

    # FBP approach
    elif ellipse_fitting.lower() == 'fbp':
        # Uses HROS, FROS, BROS as fn of WS, fuel type, etc
        b, a, c = ellipse_fbp(forward=hros_mmin,
                              flank=fros_mmin,
                              back=bros_mmin, 
                              verbose=verbose) 

    # LB provided approach
    elif lb and ellipse_fitting.lower() in ['tester', 'lb']:
        # Uses HROS and length-to-breadth ratio to calculate the ellipse
        a, b, c = ellipse_definition_tester(ros_mmin=hros_mmin, 
                                            lb=lb, 
                                            verbose=verbose)

    # Non elliptical by default (i.e., simple circle)
    else:
        a, b, c = 1., 1., 0.
        
    # Fperiods update (by default 1 per minute): e.g., fperiod_length = 2 min -> fperiods = 30 = half hour 
    if fperiod_length != 1:
        fperiods = int(fperiods / float(fperiod_length))
        print('New fperiods (due to length):', fperiods)
        
    # For fire periods
    for t in tqdm.tqdm(range(fperiods)):
        # Ignition 
        if t == 0:
            onfire.append(n_ignition)
            n_colors[n_ignition] = 'blue'

        # For on-fires
        for n_fire in onfire:
            # Get potential neighbors
            nbs = [_ for _ in G.neighbors(n_fire) if _ not in onfire]
            
            # Initialize own f_progress
            if fprogress[n_fire] == 0:
                fprogress[n_fire] = {_:0 for _ in nbs}
            

            # Update fireprogress in potential neighbors [Parallel if needed here]
            for nb in nbs:
                angle = node2angle(n_fire, nb, factor=nxm)
                if print_ros_values:
                    print('angle:', angle)
                angle, ros = angle2ros_fn(angle, thetafire, a, b, centered_at_cell=centered_at_cell)
                if print_ros_values:
                    print('\tros:', ros)
                fprogress[nb] += np.round(ros * fperiod_length, 5)
                # fprogress[n_fire][nb] += np.round(ros * fperiod_length, 5)
                if fprogress[nb] >= G.get_edge_data(n_fire, nb)['weight']:
                # if fprogress[n_fire][nb] >= G.get_edge_data(n_fire, nb)['weight']:
                    onfire.append(nb)
                    n_colors[nb] = 'red'
                    
    # Node Positions
    n_pos = {}
    for old_n in labels.keys():
        n_pos[labels[old_n]] = old_n
        
    # Create a plot
    if verbose:
        plt.figure(1,figsize=(6,6)) 
        nx.draw(G, 
                n_pos,
                font_size=8, 
                node_color=list(n_colors.values()), 
                node_size=50,
                node_shape="s",) 
        
    # Return
    return a, b, c, fprogress, onfire, n_colors, n_pos

# Get grid wrapper
def get_grid(nxm, onfire):
    # Create placeholder
    grid = np.zeros((nxm, nxm))
    
    # Fill with ones
    g_flatten = grid.flatten()
    g_flatten[onfire] = 1
    grid = g_flatten.reshape((nxm, nxm))
    
    # Rotate with the corresponding angle (90 for simplicity of the example)
    grid = np.rot90(grid)

    # Plot
    # plt.figure(1,figsize=(6,6)) 
    # sns.heatmap(grid, xticklabels=10, yticklabels=10, cbar=False)
    
    # Return binary matrix
    return grid.astype(int)