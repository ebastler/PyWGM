import sys
import os
import numpy as np
from scipy.constants import speed_of_light,pi
import re
sys.path.append(os.path.join(os.getenv('JCMROOT'), 'ThirdPartySupport', 'Python'))
import jcmwave

def runJCM(jcmMaxThreads, simFEMdeg, simPrecision, simFieldRes, modeNumber, iterGuess, iternPhi, substrateHeight, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, airHeight, absWallH, absWallW, iterPillarRadius, iterSubstrateRadius, topDBRpairs, botDBRpairs, unEtchedBotDBRpairs, permSubstrate, permBotDBR1, permBotDBR2, permCavity, permTopDBR1, permTopDBR2, permCap, permAbsWall, apertureEnable, permApertureOx, permApertureNonOx, apertureRadius, aperturePosition, apertureThickness):

    # Set priority of aperture layers to 4 and 5, respectively, when enabled (to override bot DBR)
    # Set it to 1 if not, so bot DBR will override
    if apertureEnable == True:
        apertureEnableOx = 4
        apertureEnableNonOx = 5
    else:
        apertureEnableOx = 1
        apertureEnableNonOx = 1

    # Set sim variables to pass to JCM
    keys = {
        'fieldRes'            : simFieldRes,
        'modeNumber'          : modeNumber,
        'guess'               : iterGuess,
        'Nphi'                : iternPhi,
        'pillarR'             : iterPillarRadius,
        'substrateR'          : iterSubstrateRadius,
        'substrateH'          : substrateHeight,
        'permSubstrate'       : permSubstrate,
        'cavityHeight'        : cavityHeight,
        'permCavity'          : permCavity,
        'capHeight'           : capHeight,
        'permCap'             : permCap,
        'botDBR1Height'       : botDBR1Height,
        'botDBR2Height'       : botDBR2Height,
        'topDBR1Height'       : topDBR1Height,
        'topDBR2Height'       : topDBR2Height,
        'permBotDBR1'         : permBotDBR1,
        'permBotDBR2'         : permBotDBR2,
        'permTopDBR1'         : permTopDBR1,
        'permTopDBR2'         : permTopDBR2,
        'topDBRpairs'         : topDBRpairs,
        'botDBRpairs'         : botDBRpairs,
        'airHeight'           : airHeight,
        'absWallW'            : absWallW,
        'absWallH'            : absWallH,
        'unEtchedBotDBRpairs' : unEtchedBotDBRpairs,
        'FEMD'                : simFEMdeg,
        'simPrec'             : simPrecision,
        'permAbsWallIm'       : permAbsWall,
        'apertureRadius'      : apertureRadius,
        'aperturePosition'    : aperturePosition,
        'apertureThickness'   : apertureThickness,
        'permApertureOx'      : permApertureOx,
        'permApertureNonOx'   : permApertureNonOx,
        'apertureEnableOx'    : apertureEnableOx,
        'apertureEnableNonOx' : apertureEnableNonOx,
    }
    # Write additional info to the log file
    with open("jcm_call_latest.log", "a") as file_object:
        file_object.write("\nnPhi: {}, pillarRadius: {}, guess: {}, wavelengthGuess: {}".format(iternPhi, iterPillarRadius, iterGuess, 2 * speed_of_light * pi / (iterGuess*1e-9)))
    # If jcmMaxThreads is != 0, force the set number of threads to be used
    if jcmMaxThreads > 0:
        jcmwave.set_num_threads(jcmMaxThreads)
    #
    results = jcmwave.solve('project.jcmpt', keys, logfile=open("jcm_latest.log", 'w'))
    #
    # Read complex Eigenmodes from jcm result return values, compute wavelength and Q factor
    eigenmodeRe = results[0]['eigenvalues']['eigenmode'].real
    eigenmodeIm = results[0]['eigenvalues']['eigenmode'].imag
    wavelength = 2 * speed_of_light * pi / eigenmodeRe
    qFactor = eigenmodeRe / (-2 * eigenmodeIm)
    #
    # Read the integrated field density difference between cavity and rest of the structure
    densityDiff = []
    for i in range(0, modeNumber):
        densityDiff += [results[1]['ElectricFieldEnergy'][i][6] - results[1]['ElectricFieldEnergy'][i][1] - results[1]['ElectricFieldEnergy'][i][2] - results[1]['ElectricFieldEnergy'][i][3] - results[1]['ElectricFieldEnergy'][i][4] - results[1]['ElectricFieldEnergy'][i][5]]
    #
    # Read the integrated field density in both air and structure
    energyIntAir = []
    energyIntStructure = []
    for i in range(0, modeNumber):
        energyIntAir += [results[1]['ElectricFieldEnergy'][i][0].real]
        energyIntStructure += [results[1]['ElectricFieldEnergy'][i][1].real + results[1]['ElectricFieldEnergy'][i][2].real + results[1]['ElectricFieldEnergy'][i][3].real + results[1]['ElectricFieldEnergy'][i][4].real + results[1]['ElectricFieldEnergy'][i][5].real + results[1]['ElectricFieldEnergy'][i][6].real + results[1]['ElectricFieldEnergy'][i][7].real ]
    return wavelength*1e9, qFactor, densityDiff, results[2]['field'], results[3]['field'], results[2]['X']*1e6, results[2]['Y']*1e6, results[3]['Z']*1e6, eigenmodeRe, energyIntStructure, energyIntAir

def evalJCMfields(fieldXY, fieldXZ, posX, posZ, fieldRes, wavelength, qFactor, i, permCavity):
    # Get field intensity arrays in XY (vertical cross-section) from results)
    intensityXY_X = [ [fieldXY[x][v][0].real if posX[x][v].real > 0 else -fieldXY[x][v][0].real for v in range(0, fieldRes) ] for x in range(0, fieldRes)]
    intensityXY_Y = [ [v[1].real for v in x] for x in fieldXY ]
    #
    # Check if field is polarized in X (radial) or Y (vertical) direction, use the stronger component for visualization
    intensityXY_dir = 'X' if np.amax(np.absolute(intensityXY_X)) > np.amax(np.absolute(intensityXY_Y)) else 'Y'
    intensityXY = intensityXY_X/np.amax(np.absolute(intensityXY_X)) if intensityXY_dir == 'X' else intensityXY_Y/np.amax(np.absolute(intensityXY_Y))
    #
    # Get field intensity arrays in XZ (horizontal cross-section) from results)
    # XZ_XZ field is radially polarized - this function checks for each datapoint 
    # whether the field points towards or from the origin, and assigns the length of the vector the appropriate sign
    intensityXZ_XZ = [ [ np.sqrt(np.square(fieldXZ[x][v][0].real) + np.square(fieldXZ[x][v][2].real)) if np.sqrt(np.square(fieldXZ[x][v][0].real*1e-20 + posX[x][v].real) + np.square(fieldXZ[x][v][2].real*1e-20 + posZ[x][v].real)) > np.sqrt(np.square(posX[x][v].real) + np.square(posZ[x][v].real)) else -np.sqrt(np.square(fieldXZ[x][v][0].real) + np.square(fieldXZ[x][v][2].real)) for v in range(0, fieldRes) ] for x in range(0, fieldRes) ]
    intensityXZ_Y  = [ [v[1].real for v in x] for x in fieldXZ ]
    #
    # Check if field is polarized in XZ (radial) or Y (vertical) direction, use the stronger component for visualization
    intensityXZ_dir = 'radial' if np.amax(np.absolute(intensityXZ_XZ)) > np.amax(np.absolute(intensityXZ_Y)) else 'Y'
    intensityXZ = intensityXZ_XZ/np.amax(np.absolute(intensityXZ_XZ)) if intensityXZ_dir == 'radial' else intensityXZ_Y/np.amax(np.absolute(intensityXZ_Y))
    #
    # Compute mode volume and purcell factor, script taken from another matlab script at our uni.
    # See https://thesis.library.caltech.edu/2487/5/KippenbergThesis.pdf page 24 for mode volume
    # See https://academic.oup.com/book/36428?login=false chapter 6.1 for Purcell
    ee = jcmwave.loadtable('project_results/electric_energy.jcm')
    eed = jcmwave.loadcartesianfields('project_results/electric_energy_density.jcm')
    # Find the (abs) max of the electric field
    meed = np.amax(abs(eed['field'][i]))
    # Find the X position (offset from pillar center) of the field strength max, take the abs since it is rotationally symmetric anyway and the sign is purely computation dependent
    xPosMax = abs(eed['X'][(np.where(abs(eed['field'][i]) == meed))[0]])[0]*1e9
    # Compute mode volume, no clue how this works
    modeVolume = sum((ee['ElectricFieldEnergy'][i]).real)/meed  # m^3
    modeVolumeInWL = modeVolume/wavelength**3
    # Compute purcell factor
    purcellFactor = (3/4/pi**2)*((wavelength/np.sqrt(np.real(complex(permCavity.replace('i', 'j')))))**3)*(qFactor/modeVolume)
    #
    betaFactor = purcellFactor/(purcellFactor+1)
    return intensityXZ_dir, intensityXY, intensityXZ, xPosMax, modeVolume*1e27, purcellFactor, betaFactor
