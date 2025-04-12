# coding: utf-8
__version__ = "1.0"
__author__ = "Cristobal Pais"
__collaborator__ = "Minho Kim"

# Importations
import time
import numpy as np
import numpy.random as npr
import pandas as pd
import os
import nlopt
from argparse import ArgumentParser
from Cell2Fire_DFO import Cell2Fire_Norm
from Cell2Fire_DFO import BBOFactors2CSV as BBOout

'''
Counter Global var
'''
def init():
    global counter
    counter = 0

'''
Pre-process fuels
'''
def preProcessFuels(InstancesPath, NFTypes=[-9999,0,91,92,93,98,99]):
    # Containers
    Forests = {}
    FUnique = set()

    # Instances
    for i in InstancesPath:
        print(i)
        FPath = os.path.join(i, "Forest.asc")
        Forests[i] = np.loadtxt(FPath, delimiter=" ", skiprows=6).astype(np.int32)
        print(np.unique(Forests[i]))
        FUnique |= set(np.unique(Forests[i]))

    for i in NFTypes:
        if i in FUnique:
            FUnique.remove(i)
            
    print("Unique fuels:", FUnique) if verbose else None
    print("Total number of unique fuels:", len(FUnique)) if verbose else None
    return FUnique

'''
BBO Factors file generator
Creates .csv file to be read from C++ with factors
'''
def BBOFactors2CSV(OutPath, Factors=None, nfactors=4, verbose=False):
    # Create BBOFuel.csv file (or replace an old one)
    DFcomp = {}
    FactorList = np.ones(nfactors)
    row = 0
    if verbose:
        print("Unique Fuels:", FUnique)
    for ftype in FUnique:
        if Factors is not None:
            for aux in range(0, nfactors):
                FactorList[aux] = Factors[aux + 4*row]
            row += 1
        DFcomp[ftype] = FactorList.copy()
    DF = pd.DataFrame(DFcomp)
    DF = DF.transpose()
    DF.index.name = "FType"
    DF = DF.rename(columns={0: "HFactor", 1:"FFactor", 2:"BFactor", 3:"EFactor"})
    DF.to_csv(os.path.join(OutPath, "BBOFuels.csv"), header=True,  sep=",")
    
    # Return to fill vector X
    return DF.values.flatten()

    

# Main
if __name__ == "__main__":
    # Read args
    parser = ArgumentParser()
    parser.add_argument("--PathC",
                        help="Path to the main exe file (Cell2Fire Sim)",
                        dest="PathC",
                        type=str,
                        default=None)
    parser.add_argument("--OutFolder",
                        help="Path to the output folder (Cell2Fire Sim)",
                        dest="OutFolder",
                        type=str,
                        default=None)
    parser.add_argument("--INSTANCE",
                        help="Instance name for reference fire simulation output",
                        dest="INSTANCE",
                        type=str,
                        default="Sim")                        
    parser.add_argument("--PathScars",
                        help="Path to the original scar folder (Farsite/Real fire)",
                        dest="PathScars",
                        type=str,
                        default=None)
    parser.add_argument("--input-instance-folder",
                        help="Path to the instance folder (Cell2Fire Sim)",
                        dest="PathInstance",
                        type=str,
                        default=None)
    parser.add_argument("--fobjective",
                        help="Objective function for the optimization problem",
                        dest="fobjective",
                        type=str,
                        default="norm")
    parser.add_argument("--ROS_Threshold",
                        help="SROS critical threshold (m/s) default 0.1.",
                        dest="ROS_Threshold",
                        type=float,
                        default=0.1)      
    parser.add_argument("--HFI_Threshold",
                        help="Head Fire Intensity critical threshold (m/s) default 0.1.",
                        dest="HFI_Threshold",
                        type=float,
                        default=0.1)      
    parser.add_argument("--norm",
                        help="Norm for the fire scar function (Cell2Fire Sim)",
                        dest="normF",
                        type=str,
                        default="fro")
    parser.add_argument('--muWeights', nargs='+', 
                        help="Hour weights (per hour fire scar weight)",
                        type=float, 
                        dest="muWeights",
                        default=np.ones(10))
    parser.add_argument('--toTune', nargs='+', 
                        help="Components to tune, combination of: SROS, CriticalSROS, CROS, CriticalCROS, RTime, EllipticalROS",
                        type=str, 
                        dest="toTune",
                        default="")
    parser.add_argument("--UB",
                        help="Variables upper-bound",
                        dest="UB",
                        type=float,
                        default=10.)
    parser.add_argument("--LB",
                        help="Variables lower-bound",
                        dest="LB",
                        type=float,
                        default=0.)
    parser.add_argument("--absTol",
                        help="Absolute tolerance (stopping criterion)",
                        dest="absTol",
                        type=float,
                        default=1e-10)  
    parser.add_argument("--algorithm",
                        help="Select the algorithm to use (default BOBYQA): COBYLA, BOBYQA, NEWUOA_BOUND, LN_NEWUOA, LN_PRAXIS, LN_NELDERMEAD, LD_MMA, LD_SLSQP, LD_LBFGS, AUGLAG",
                        dest="alg",
                        type=str,
                        default='BOBYQA')
    parser.add_argument("--verbose",
                        help="Verbosity level",
                        dest="verbose",
                        action='store_true',
                        default=False)
    parser.add_argument("--cros",
                        help="Allow CROS during the simulation",
                        dest="cros",
                        action='store_true',
                        default=False)
    parser.add_argument("--randomInit",
                        help="Random initialization of the factors",
                        dest="randomInit",
                        action='store_true',
                        default=False)
    parser.add_argument("--seed",
                        help="Random seed",
                        dest="seed",
                        type=int,
                        default=123)
    parser.add_argument("--randomConstant",
                        help="Constant for random initialization",
                        dest="rndconstant",
                        type=int,
                        default=5)
    parser.add_argument("--saveFreq",
                        help="Save grid and weights of the current iteration (-1 never)",
                        dest="saveFreq",
                        type=int,
                        default=-1)
    
    # Parse arguments
    args = parser.parse_args()
    
    verbose = args.verbose
    INSTANCE = args.INSTANCE if args.INSTANCE else os.path.basename(args.PathInstance)

    # Global counter
    init()
    
    # Flags (inside dictionary)
    TuningParams = {}
    LengthParams = np.zeros(2).astype(np.int32)
    
    # Fuel types
    FUnique = preProcessFuels([args.PathInstance], NFTypes=[-9999,0,91,92,93,98,99])
    
    # Initialization strategies
    if args.randomInit is False:
        # Non-Random initialization
        # Initialize BBOFuels.csv if not existent or read an existing one (starting point)
        BBOFile = os.path.join(args.PathInstance, "BBOFuels.csv")
        if os.path.isfile(BBOFile) is False:
            toFillX = BBOFactors2CSV(args.PathInstance, Factors=None, nfactors=4, verbose=False)
        else:
            toFillX = pd.read_csv(os.path.join(args.PathInstance, "BBOFuels.csv"), sep=",", index_col="FType").values.flatten()
    else:
        # Random factors
        toFillX = args.rndconstant * npr.random(len(FUnique) * 4)
    
    # Filter adjustement parameters
    if len(args.toTune) > 0:
        # Sort toTune
        sortedTune = []
        if "CriticalSROS" in args.toTune:
            sortedTune.append("CriticalSROS")
        if "EllipticalROS" in args.toTune:
            sortedTune.append("EllipticalROS")
        if "CriticalCROS" in args.toTune:
            sortedTune.append("CriticalCROS")
        if "CriticalActiveCROS" in args.toTune:
            sortedTune.append("CriticalActiveCROS")
        if "CROS_FM10" in args.toTune:
            sortedTune.append("CROS_FM10")
        if "CROS_CBD" in args.toTune:
            sortedTune.append("CROS_CBD")
        if "CROS_CCF" in args.toTune:
            sortedTune.append("CROS_CCF")        
        
        # Initialize x
        x_0 = np.array([])
        for i in sortedTune:
            if i == "CriticalSROS":
                x_0 = np.concatenate((x_0, [args.ROS_Threshold]))
                TuningParams[i] = True
                LengthParams[0] = 1
                LengthParams[1] = LengthParams[0].copy()
            if i == "EllipticalROS":
                x_0 = np.concatenate((x_0, np.asarray(toFillX)))
                TuningParams[i] = True
                LengthParams[1] += 4 * len(FUnique)             
    
    # Information
    print("TuningParams:", TuningParams) if verbose else None
    print("LengthParams:", LengthParams) if verbose else None
         
    # NLOPT package
    # Info
    print("Total parameters: 1 Critical SROS, 4 EllipticalROS *", len(FUnique), "FuelTypes") if verbose else None
    print("Tuning:", [k for k in TuningParams.keys()]) if verbose else None
    print("X0:", x_0, "\nN:", len(x_0)) if verbose else None
    time.sleep(2)

    # BOBYQA
    print("\n\n********* Function Cell2Fire  ********") if verbose else None
    alg = "BOBYQA"
    print("Algorithm:", alg) if verbose else None

    # Dimension
    n = len(x_0)                           
    
    '''
    NLOPT_LN_COBYLA
    NLOPT_LN_BOBYQA
    NLOPT_LN_NEWUOA_BOUND  |  NLOPT_LN_NEWUOA
    NLOPT_LN_PRAXIS
    NLOPT_LN_NELDERMEAD
    NLOPT_LD_MMA
    NLOPT_LD_SLSQP
    NLOPT_LD_LBFGS
    NLOPT_AUGLAG
    '''
    
    # Create the opt object
    opt = nlopt.opt(nlopt.LN_BOBYQA, n)

    # Bounds
    opt.set_lower_bounds(np.full(n, args.LB))
    opt.set_upper_bounds(np.full(n, args.UB))

    # Checks
    if np.any(x_0 < args.LB) or np.any(x_0 > args.UB):
        print("WARNING: Initial point outside bounds!")
        # Clip initial point to bounds
        x_0 = np.clip(x_0, args.LB, args.UB)
        # Add before optimizer setup
        print(f"x_0 dimension: {len(x_0)}")
        print(f"Optimizer dimension: {n}")
        print(f"Bounds arrays dimensions: {len(np.full(n, args.LB))}, {len(np.full(n, args.UB))}")        
    # Make sure LB is strictly less than UB
    if args.LB >= args.UB:
        print("ERROR: Lower bound must be less than upper bound")
        args.LB = args.UB - 0.1  # Adjust as needed
        # Right before the optimize call, add:
        print("Algorithm:", opt.get_algorithm_name())
        print("Dimension:", opt.get_dimension())
        print("Lower bounds:", opt.get_lower_bounds())
        print("Upper bounds:", opt.get_upper_bounds())

    # Objective
    opt.set_min_objective(lambda x, grad: Cell2Fire_Norm(x, grad, args.OutFolder, args.PathC, 
                                                         args.PathScars, args.PathInstance,
                                                         args.muWeights, FUnique, INSTANCE, 
                                                         args.normF, TuningParams, LengthParams,
                                                         args.ROS_Threshold, args.HFI_Threshold, args.cros,
                                                         args.fobjective, args.saveFreq, args.verbose))
    
    # Tolerance
    opt.set_xtol_abs(args.absTol)
    #opt.set_xtol_rel(1e-26)
    
    # Time and optimization
    t1= time.time()
    x = opt.optimize(x_0)
    t2 = time.time()
    
    # Minimum fobj
    minf = opt.last_optimum_value()

    # Final results
    print ("\n")
    print ("---------- Results from BOBYQA -----------------")
    print("Optimum at:\t", x)
    _ = BBOout(args.PathInstance, FUnique, x.copy())
    print("Minimum obj. value:\t", minf)
    print("Number of evaluations:\t", opt.get_numevals())
    print("Total Runtime [s]:\t", t2-t1)
