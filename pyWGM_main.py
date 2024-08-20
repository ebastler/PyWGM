# Moritz Plattner, 2021
# Python script to simulate various properties of WGM modes in micropillars
# This script is based on JCMwave and needs a valid install + licence in order to work!
# It calls multiple functions from other files, please do not rename or move python files in this directory.

# Import necessary libs, not all needed in this file but listed anyway
# import sys
# import os
# import csv
# import re
import math
# sys.path.append(os.path.join(os.getenv('JCMROOT'), 'ThirdPartySupport', 'Python'))
# import jcmwave
# import matplotlib.pyplot as plt
# from scipy.constants import speed_of_light,pi
import time


# Choose which sims to run (True = yes, False = no)
# Warning: Running all sims can take a long time, depending on parameter choice, and is usually not reasonable.
runDiaScan      = True                                      # Vary the pillar diameter
runApertureScan = False                                     # Simulate the same structure, once with various aperture sizes and without aperture
runTopDbrScan   = False                                     # Vary the number of top DBR mirror pairs
runDbrEtchScan  = False                                     # Vary the amount of etched bottom DBR pairs
runRefIndexScan = False                                     # Vary the refractive index of the cavity

# Set this line to the number of threads your PC/VM has. Needed on some machines to use 100% of the CPU. Set to 0 to let JCM auto-detect.
jcmMaxThreads = 0

# Set material parameters
refIndexGaAs = 3.5480
refIndexAl90Ga10As = 3.0150
refIndexAl45Ga55As = 3.3041
refIndexAl20Ga80As = 3.4225
refIndexAlAs = 2.9434
refIndexAL2O3 = 1.7565

refIndexSiN = 1.886
refIndexSiO = 1.9108
# Compute relative permittivity from ref index, or enter manually
permGaAs = refIndexGaAs**2
permAl20Ga80As = refIndexAl20Ga80As**2
permAl45Ga55As = refIndexAl45Ga55As**2
permAl90Ga10As = refIndexAl90Ga10As**2
permAlAs = refIndexAlAs**2
permAl2O3 = refIndexAL2O3**2
permSiN = refIndexSiN**2
permSiO = refIndexSiO**2

# Set permittivities for different parts of the structure
permBotDBR1 = permAl90Ga10As                                # First layer above substrate
permBotDBR2 = permAl20Ga80As                                # Second layer above substrate
permTopDBR1 = permAl90Ga10As                                # First layer above cavity
permTopDBR2 = permAl20Ga80As                                # Second layer above cavity
permSubstrate = permGaAs                                    # Permittivity of the substrate
permCap = permGaAs                                          # Permittivity of the capping layer
permApertureOx = permAl2O3                                  # Permittivity of the oxidized part of the AlAs layer
permApertureNonOx = permAlAs                                # Permittivity of the non oxidized part of the AlAs layer

# Old values used for previous cals - left here for comparison
permAbsWall = '{}+.7i'.format(permGaAs)                     # Permittivity of the absorbing wall - should be complex to add absorption.
permCavity = '{}+0.0005i'.format(permGaAs)                  # Permittivity of the cavity - should be complex to simulate absorption of QDs and surface states

# Set structure parameters [nm]
minPillarRadius = 1000                                      # lowest radius to be simulated in nm
maxPillarRadius = 2500                                      # max radius to be simulated in nm

substrateHeight = 300                                       # height of the substrate
cavityHeight = 275.4                                        # cavity thickness
capHeight = 53                                              # height of the GaAs top capping layer

botDBR1Height = 83.9                                        # First layer above substrate
botDBR2Height = 73.7                                        # Second layer above substrate
topDBR1Height = 83.9                                        # First layer above cavity
topDBR2Height = 73.7                                        # Second layer above cavity

topDBRpairs = 12                                            # number of top DBR pairs
botDBRpairs = 32                                            # number of bottom DBR pairs

absWallW = 2                                                # thickness of the absorbing wall on the pillar outside to simulate roughness.

unEtchedBotDBRpairs = math.ceil(botDBRpairs/2)              # number of un-etched bottom DBR pairs (default: half, ceiling function for odd numbers)

apertureEnable = True                                       # enable or disable aperture for site controlled growth
minApertureSize = 180                                       # Size of the oxidized ring - 200 = 200nm smaller than pillar radius. Min value for runApertureScan. 
maxApertureSize = 250                                       # Only used for for runApertureScan
apertureThickness = 30                                      # Thickness of the AlAs layer that will be oxidized
apertureSpacing = 40                                        # Distance between aperture (top) and cavity (bottom)

# Set variation limits for simulations
# dbrEtchScan
minUnEtchedBotDBRpairs = 0                                  # lowest number of un-etched bottom DBR pairs (default: 0)
maxUnEtchedBotDBRpairs = botDBRpairs                        # highest number of un-etched bottom DBR pairs (default: all)

# dbrDetScan
layerheightVar = 0.6                                        # vary DBR thickness between (1-x) and (1+x) times the target thickness
layerheightVarStep = 0.1                                    # step size for layerheight variation

# diaScan
diaScanSteps = 250                                          # steps between radius in nm

# topDbrScan
minTopDBRpairs = 12
maxTopDBRpairs = 16

# apertureScan
apertureSizeStepping = 50

# runRefIndexScan
minRefindexCavity = 3.4
maxRefindexCavity = 4.2
refindexCavityStep = 0.2

# Simulated wavelength ranges
startWL = 930                                               # Shortest wavelentgh to simulate
stopWL = 990                                                # Longest wavelentgh to simulate
targetWL = 960                                              # WL for which you want the FSR calculated - first mode above and below is used for FSR.

# Set sim parameters
simFEMdeg = 3                                               # FEM degreee of the simulation
simPrecision = 5e-4                                         # Precision of the simulation
simFieldRes = 512                                           # number of datapoints to export per field and axis (export res for field plot visualisation only)
modeNumber = 2                                              # number of modes for JCM to compute. 2 is fine if your guess is good, increase if you do not find 2 modes per step


# Start time measurment
startTime = time.time()

if runDiaScan:
    import pyWGM_scandia
    pyWGM_scandia.rundiascan(substrateHeight, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, minPillarRadius, maxPillarRadius, diaScanSteps, absWallW, topDBRpairs, botDBRpairs, unEtchedBotDBRpairs, minApertureSize, apertureThickness, apertureSpacing, apertureEnable, permApertureOx, permApertureNonOx, permSubstrate, permBotDBR1, permBotDBR2, permCavity, permTopDBR1, permTopDBR2, permCap, permAbsWall, simFEMdeg, simPrecision, simFieldRes, modeNumber, jcmMaxThreads, startWL, stopWL, targetWL)

if runApertureScan:
    import pyWGM_scanaperture
    pyWGM_scanaperture.runaperturescan(substrateHeight, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, minPillarRadius, absWallW, topDBRpairs, botDBRpairs, unEtchedBotDBRpairs, minApertureSize, maxApertureSize, apertureSizeStepping, apertureThickness, apertureSpacing, apertureEnable, permApertureOx, permApertureNonOx, permSubstrate, permBotDBR1, permBotDBR2, permCavity, permTopDBR1, permTopDBR2, permCap, permAbsWall, simFEMdeg, simPrecision, simFieldRes, modeNumber, jcmMaxThreads, startWL, stopWL, targetWL)

if runTopDbrScan:
    import pyWGM_scantopdbr
    pyWGM_scantopdbr.runtopdbrscan(substrateHeight, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, minPillarRadius, absWallW, minTopDBRpairs, maxTopDBRpairs, botDBRpairs, unEtchedBotDBRpairs, minApertureSize, apertureThickness, apertureSpacing, apertureEnable, permApertureOx, permApertureNonOx, permSubstrate, permBotDBR1, permBotDBR2, permCavity, permTopDBR1, permTopDBR2, permCap, permAbsWall, simFEMdeg, simPrecision, simFieldRes, modeNumber, jcmMaxThreads, startWL, stopWL, targetWL)

if runDbrEtchScan:
    import pyWGM_scandbretch
    pyWGM_scandbretch.rundbretchscan(substrateHeight, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, minPillarRadius, absWallW, topDBRpairs, botDBRpairs, unEtchedBotDBRpairs, minUnEtchedBotDBRpairs, maxUnEtchedBotDBRpairs, minApertureSize, apertureThickness, apertureSpacing, apertureEnable, permApertureOx, permApertureNonOx, permSubstrate, permBotDBR1, permBotDBR2, permCavity, permTopDBR1, permTopDBR2, permCap, permAbsWall, simFEMdeg, simPrecision, simFieldRes, modeNumber, jcmMaxThreads, startWL, stopWL, targetWL)

if runRefIndexScan:
    import pyWGM_scanrefindex
    pyWGM_scanrefindex.runrefindexscan(substrateHeight, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, minPillarRadius, absWallW, topDBRpairs, botDBRpairs, unEtchedBotDBRpairs, minApertureSize, apertureThickness, apertureSpacing, apertureEnable, permApertureOx, permApertureNonOx, permSubstrate, permBotDBR1, permBotDBR2, permCavity, minRefindexCavity, maxRefindexCavity, refindexCavityStep, permTopDBR1, permTopDBR2, permCap, permAbsWall, simFEMdeg, simPrecision, simFieldRes, modeNumber, jcmMaxThreads, startWL, stopWL, targetWL)


print("The sim run took {} seconds".format(time.time() - startTime))