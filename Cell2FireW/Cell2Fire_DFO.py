# coding: utf-8
__version__ = "1.0"
__author__ = "Cristobal Pais"
__collaborator__ = "Minho Kim"

# Importations
import numpy as np
import pandas as pd
import os
import shutil
import time
import subprocess
import rasterio

# Deal with Paths
import os.path

# Scipy for linear algebra
from scipy.sparse import *
from scipy import sparse
from scipy.sparse.linalg import norm
from skimage.metrics import structural_similarity, mean_squared_error

# No warnings
import warnings
warnings.simplefilter("error", RuntimeWarning)

# Global counter
global counter
counter = 0

'''
BBO Factors file generator
Creates .csv file to be read from C++ with factors
'''
def BBOFactors2CSV(OutPath, FUnique, Factors, nfactors=4, verbose=False, fname="BBOFuels.csv"):
    # Create BBOFuel.csv file (or replace an old one)
    DFcomp = {}
    FactorList = np.ones(nfactors)
    row = 0
    # print("FUnique:", FUnique)
    for ftype in FUnique:
        for aux in range(0, nfactors):
            FactorList[aux] = Factors[aux + 4*row]
        row += 1
        DFcomp[ftype] = np.round(FactorList, 3).copy()
    DF = pd.DataFrame(DFcomp)
    DF = DF.transpose()
    DF.index.name = "FType"
    DF = DF.rename(columns={0: "HFactor", 1:"FFactor", 2:"BFactor", 3:"EFactor"})
    DF.to_csv(os.path.join(OutPath, fname), header=True,  sep=",")
    
'''
Function for BOBYQUA
Description:
Call to the simulator norm difference function
'''
def Cell2Fire_Norm(x, grad, OutFolder, PathC, PScars, PInstance, weights, FUnique, INSTANCE, ntype='fro',
                   TuningParams={}, LengthParams={}, ROS_Threshold=0.0, HFI_Threshold=0.0, cros=False, 
                   fobjective="norm", saveFreq=-1, verbose=False):
    # Check global var
    global counter 
    counter += 1

    # X point
    if grad.size > 0:
        pass
    
    x = np.abs(x).flatten()
    print("****************************************************************************************") if verbose else None
    print("Iteration:", counter)
    print("Objective function:", fobjective) if verbose else None
    print("\nInitial x:", x) if verbose else None  
       
    #Cell2Fire_Norm function
    PathOut = OutFolder
    PathScars = PScars
    PathCSim = PathC
    PathInstance = PInstance
    
    # Objective function
    Fobj = 0
    
    # Command line 
    # Parse args for calling C++ via subprocess        
    execArray=[os.path.join(os.getcwd(),'Cell2Fire'), 
               '--input-instance-folder', PathInstance,
               '--output-folder', OutFolder,
               '--ignitions',
               '--final-grid',
               '--nsims', str(1),
               '--sim-years', str(1),
               '--weather', 'rows',
               '--nweathers', str(1),
               '--sim', 'US',
            #    '--Fire-Period-Length', str(1.0),
               '--ROS-Threshold' if ("CriticalSROS" in TuningParams.keys()) else '',   
               str(x[0]) if ("CriticalSROS" in TuningParams.keys()) else '',
               '--HFI-Threshold', str(HFI_Threshold),
               '--CROS-Threshold' if ("CriticalCROS" in TuningParams.keys()) else '',   
               str(x[1]) if ("CriticalCROS" in TuningParams.keys()) else '',
               '--CROS-Active-Threshold' if ("CriticalActiveCROS" in TuningParams.keys()) else '',   
               str(x[2]) if ("CriticalActiveCROS" in TuningParams.keys()) else '',
               '--CROS-CBD' if ("CROS_CBD" in TuningParams.keys()) else '',   
               str(x[2]) if ("CROS_CBD" in TuningParams.keys()) else '',
               '--CROS-CCF' if ("CROS_CCF" in TuningParams.keys()) else '',   
               str(x[2]) if ("CROS_CCF" in TuningParams.keys()) else '',
               '--CROS-FM10' if ("CROS_FM10" in TuningParams.keys()) else '',   
               str(x[2]) if ("CROS_FM10" in TuningParams.keys()) else '',
               '--cros' if (cros is True) else '',
               '--bbo',
            #    '--verbose'
              ]
    
    if verbose:
        print("ExecArray:", execArray)

    # Information
    if "CriticalSROS" in TuningParams.keys():
        if "EllipticalROS" in TuningParams.keys():      
            BBOFactors2CSV(PInstance, FUnique, x.copy())
            
    else:
        BBOFactors2CSV(PInstance, FUnique, x.copy())
            
    # Output log
    if OutFolder is not None:
        if os.path.isdir(OutFolder) is False:
            os.makedirs(OutFolder)
    LogName = os.path.join(OutFolder, "LogFile.txt")

    # Perform the call
    with open(LogName, 'w') as output:
        proc = subprocess.Popen(execArray, stdout=output)
        proc.communicate()
    return_code = proc.wait()
    if (return_code != 0):
        raise RuntimeError(f'C++ returned {return_code}.\nTry looking at {LogName}.') 
        

    # Perform the call
    try:
        with open(LogName, 'w') as output:
            proc = subprocess.Popen(execArray, stdout=output)
            proc.communicate()
        return_code = proc.wait()
        if (return_code != 0):
            raise RuntimeError(f'C++ returned {return_code}.\nTry looking at {LogName}.') 
        print("done running sim...") if verbose else None
                
        # Reading original scar
        GridPath = os.path.join(OutFolder, "Grids", "Grids1")
        if os.path.exists(GridPath):

            GridFiles = os.listdir(GridPath)
            print("Info:", GridPath, GridFiles) if verbose else None
            if len(GridFiles) > 0: 
                ForestGridM = pd.read_csv(GridPath +"/"+ GridFiles[-1], delimiter=',', header=None).values
                print("ForestGridM shape:", ForestGridM.shape) if verbose else None
                sM1 = sparse.csr_matrix(ForestGridM)
                # sM1 = ForestGridM.flatten()
                
                # Comparison grids -- Load reference output to calculate error 
                # /Users/minho/Documents/GitHub/Cell2FireML/data/farsite/sherpa_10H.tif
                filenameScar = os.path.join(PathScars, INSTANCE + ".tif")
                # RealForestGridM = pd.read_csv(filenameScar, header=None).values
                RealForestGridM = rasterio.open(filenameScar).read(1)
                RealForestGridM[RealForestGridM<=0] = 0
                RealForestGridM[RealForestGridM>0] = 1

                print("RealForestGridM shape:", RealForestGridM.shape) if verbose else None
                sM2 = sparse.csr_matrix(RealForestGridM)
                # sM2 = RealForestGridM.flatten()

                # Compute the norm
                if fobjective == "norm":
                    Fobj+= norm(sM1 - sM2, ntype)
                if fobjective == "hybrid_norm":
                    Fobj+= 100. * np.abs(np.sum(sM1) - np.sum(sM2)) / np.sum(sM1) + norm(sM1 - sM2, ntype)
                if fobjective == "mse":
                    Fobj+= mean_squared_error(RealForestGridM, ForestGridM) 
                if fobjective == "ssim":
                    Fobj+= - structural_similarity(RealForestGridM, ForestGridM, multichannel=None, data_range=1)  
                if fobjective == "largeFires":
                    Delta = RealForestGridM - ForestGridM
                    Fobj+= 0.5 * np.sum(Delta[Delta == 1]) - 1.0 * np.sum(Delta[Delta == -1])
                if fobjective == "smallFires":
                    Delta = RealForestGridM - ForestGridM
                    Fobj+= 1.0 * np.sum(Delta[Delta == 1]) - 0.5 * np.sum(Delta[Delta == -1])
                if fobjective == "absolute":
                    Fobj+= np.sum(np.abs(sM1-sM2)).astype(np.float32)
                if fobjective == "burned":
                    Fobj+= np.abs(np.sum(RealForestGridM[RealForestGridM == 1]) - np.sum(ForestGridM[ForestGridM == 1])).astype(np.float)
                
                # Save grid and BBO weights if asked
                if saveFreq > 0 and (counter % saveFreq == 0 or counter == 1):
                    print("Saving the current grid and parameters for tracking the evolution of the adjustment") if verbose else None
                    instance = PathOut.split("\\")[-1]
                    SavePath = os.path.join(OutFolder, instance + "_Evolution")

                    # Check if path exists, if not, create it
                    if os.path.exists(SavePath) is False:
                        os.mkdir(SavePath)
                    
                    fileName = os.path.join(SavePath, "FireScar_" + str(counter) + ".csv")
                    fileNameBBO = "BBOFuels_" + str(counter) + ".csv"
                    
                    # Save the grid and weights
                    np.savetxt(fileName, ForestGridM.astype(np.int32), fmt="%i", delimiter=" ")
                    BBOFactors2CSV(SavePath, FUnique, x.copy(), fname=fileNameBBO)
                    
                    # Save the grid
                    np.savetxt(fileName, ForestGridM.astype(np.int32), fmt="%i", delimiter=" ")

                # Test objective function 
                shutil.rmtree(PathOut + "/Grids") 

        else:
            print("No grids were generated from the simulation, check parameters...")
            Fobj = 99999.
        
    
    except Exception as e:
        # No solution, then, return inf
        print("Algorithm failed due to:", e)
        Fobj = 99999.
        
    
    # Delta time
    print("Fobj:", Fobj)
    time.sleep(1)
    
    # Return
    return Fobj