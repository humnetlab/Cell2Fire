# Utilities
import os, subprocess
import rasterio 
import numpy as np
import pandas as pd

import os, glob, time
import matplotlib.pyplot as plt

# Basic metrics
from numpy import linalg as LA
from sklearn.metrics import f1_score
from skimage.metrics import structural_similarity as ssim


# Print solutions from results folders
def print_grids_results_folders(folders):
    # Sort
    folders.sort()
    
    # Load and print
    for f in folders:
        grids = glob.glob(os.path.join(f, 'Grids', 'Grids1', '*.csv'))
        file = grids[0]

        # grid = np.loadtxt(file, delimiter=' ')
        grid = np.loadtxt(file, delimiter=',')

        # Plot
        plt.figure(1, figsize=(6,6)) 
        plt.imshow(grid)
        plt.title(f.split(os.path.sep)[-1])
        plt.show()
        time.sleep(0.1)
    
    # Return grid
    return grid
        
# Print solutions from Farsite
def print_grids_farsite_folders(folders, filename='', delimiter=' '):
    # Sort
    folders.sort()
    
    # Load and print
    for f in folders:
        grids = glob.glob(os.path.join(f, filename))
        file = grids[0]

        # Space delimiters
        grid = rasterio.open(file)
        grid = grid.read()
        grid = grid[0, :, :]
        
        grid[grid > 0] = 1
        grid[grid == -9999] = 0
        grid = grid.astype(int)

        # Plot
        plt.figure(1, figsize=(6,6)) 
        plt.imshow(grid, cmap='Reds')
        plt.title(f.split(os.path.sep)[-1])
        plt.show()
        time.sleep(0.1)
    
    # Return grid
    return grid
    
# Open log and print
def check_log(args):
    with open(os.path.join(args.OutFolder, "LogFile.txt"), "r") as f:
        file = f.readlines()
    print(''.join(file))
    
# Update paths from args
def update_args_paths(args, INSTANCE_PATH, OUT_PATH, instance_name, out_name):
    args['InFolder'] = INSTANCE_PATH + instance_name + '/'
    args['OutFolder'] = OUT_PATH + out_name + '/'
    
    # Convert to object
    args_obj = Dict2Class(args)
    
    # Return 
    return args_obj

# Update paths from args
def update_args_full_paths(args, INSTANCE_PATH, OUT_PATH, instance_name, out_name):
    args['InFolder'] = INSTANCE_PATH + instance_name + '/'
    args['OutFolder'] = OUT_PATH + out_name + '/'
    
    # Convert to object
    args_obj = Dict2Class(args)
    
    # Return 
    return args_obj

# Print solutions from kitral
def print_grids_kitral_folders(folders, filename='kitral_scar.txt', delimiter=' '):
    # Sort
    folders.sort()
    
    # Load and print
    for f in folders:
        grids = glob.glob(os.path.join(f, filename))
        file = grids[0]

        # Space delimiters
        grid = np.loadtxt(file, delimiter=delimiter)

        # Plot
        plt.figure(1, figsize=(6,6)) 
        plt.imshow(grid, cmap='Reds')
        plt.title(f.split(os.path.sep)[-1])
        plt.show()
        time.sleep(0.1)
    
    # Return grid
    return grid
    
# Print solutions from prometheus
def print_grids_prometheus_folders(folders, filename='PromGrid7.txt', delimiter=' '):
    # Sort
    folders.sort()
    
    # Load and print
    for f in folders:
        grids = glob.glob(os.path.join(f, filename))
        file = grids[0]

        # Space delimiters
        grid = np.loadtxt(file, delimiter=delimiter)

        # Plot
        plt.figure(1, figsize=(6,6)) 
        plt.imshow(grid, cmap='Reds')
        plt.title(f.split(os.path.sep)[-1])
        plt.show()
        time.sleep(0.1)
    
    # Return grid
    return grid

# Print solutions from Farsite
def print_grids_farsite_folders(folders, filename='', delimiter=' '):
    # Sort
    folders.sort()
    
    # Load and print
    for f in folders:
        grids = glob.glob(os.path.join(f, filename))
        file = grids[0]

        # Space delimiters
        grid = rasterio.open(file)
        grid = grid.read()
        grid = grid[0, :, :]
        
        grid[grid > 0] = 1
        grid[grid == -9999] = 0
        grid = grid.astype(int)

        # Plot
        plt.figure(1, figsize=(6,6)) 
        plt.imshow(grid, cmap='Reds')
        plt.title(f.split(os.path.sep)[-1])
        plt.show()
        time.sleep(0.1)
    
    # Return grid
    return grid

# Metrics
def calculate_metrics(cell2fire, scar, instance_name='', fspread_name='Kitral'):
    norm = LA.norm(cell2fire - scar)
    mse = ((cell2fire - scar)**2).mean()
    f1 = f1_score(scar.flatten() , cell2fire.flatten(), average='macro')
    ssim_sc = ssim(cell2fire, scar, data_range=1)
    area_cell2fire = len(cell2fire[cell2fire == 1])
    area_scar = len(scar[scar == 1])

    # Results
    df_results = pd.DataFrame({'Instance':[instance_name],
                              '$delta$ norm':[norm],
                              'MSE': [mse],
                              'F1': [f1],
                              'SSIM': [ssim_sc],
                              'AreaCell2Fire': [area_cell2fire],
                              'Area' + fspread_name: [area_scar],
                              }).round(3)
    return df_results

    
# Turns a dictionary into a class
class Dict2Class(object):
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])
    def add_entries(self, entries):
        for key, value in entries.items():
            setattr(self, key, value)  

# Pre-processing: Generate the Data.csv file for the C++ core
def generateDataC(args):
    dataName = os.path.join(args.InFolder, "Data.csv")
    if os.path.isfile(dataName) is False:
        print("Generating Data.csv File...")
    if args.region == "canada":
        import DataGeneratorC as DataGenerator
        DataGenerator.GenDataFile(args.InFolder) 
    else:
        import DataGeneratorC_Chile as DataGenerator
        DataGenerator.GenDataFile(args.InFolder) 
   
                    
# Run C++ Sim 
def run(cell2fire_path, args, opt_args=None, verbose=False):
    # Parse args for calling C++ via subprocess        
    # execArray=[os.path.join(os.getcwd(),'Cell2FireC/Cell2Fire'), 
    execArray=[os.path.join(cell2fire_path,'Cell2Fire'), 
               '--input-instance-folder', args.InFolder,
               '--output-folder', args.OutFolder if (args.OutFolder is not None) else '',
               '--ignitions' if (args.ignitions) else '',
               '--sim-years', str(args.sim_years),
               '--nsims', str(args.nsims),
               '--grids' if (args.grids) else '', 
               '--final-grid' if (args.finalGrid) else '',
               '--Fire-Period-Length', str(args.input_PeriodLen),
               '--output-messages' if (args.OutMessages) else '',
               '--weather', args.WeatherOpt,
               '--nweathers', str(args.nweathers),
               '--ROS-CV', str(args.ROS_CV),
               '--IgnitionRad', str(args.IgRadius), 
               '--seed', str(int(args.seed)),
               '--nthreads', str(int(args.nthreads)),
               '--ROS-Threshold', str(args.ROS_Threshold),
               '--HFI-Threshold', str(args.HFI_Threshold),
               '--HFactor', str(args.HFactor),
               '--BFactor', str(args.BFactor),
               '--FFactor', str(args.FFactor),
               '--EFactor', str(args.EFactor),
               '--bbo' if (args.BBO) else '',
            #    '--deactivate-offset', args.offset,
               '--HarvestPlan', args.HCells if(args.HCells is not None) else '',
               '--verbose' if (args.verbose) else '',]


    # Ellipse optimization for US (FarSite)
    if opt_args:
        execArray += [
            '--kopt', str(opt_args['kopt']),
            '--EllipticalOption', str(opt_args['EllipticalOption']),
            '--LBFormula', str(opt_args['LBFormula']),
        ]

    # Print exec array
    if verbose:
        print('ExecArray:', ' '.join(execArray))
    
    # Output log
    if args.OutFolder is not None:
        if os.path.isdir(args.OutFolder) is False:
            os.makedirs(args.OutFolder)
        LogName = os.path.join(args.OutFolder, "LogFile.txt")
    else:
        LogName = os.path.join(args.InFolder, "LogFile.txt")   

    # Perform the call
    with open(LogName, 'w') as output:
        proc = subprocess.Popen(execArray, stdout=output)
        proc.communicate()
    return_code = proc.wait()
    if (return_code != 0):
        raise RuntimeError(f'C++ returned {return_code}.\nTry looking at {LogName}.') 

    # End of the replications
    print("End of Cell2FireC execution...")

    return execArray
