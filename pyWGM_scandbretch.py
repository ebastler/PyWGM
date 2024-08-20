import numpy as np
import os
import csv
import pyWGM_plot
import pyWGM_runjcm
import pyWGM_modeestimation
from scipy.constants import speed_of_light,pi

def rundbretchscan(substrateHeight, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, minPillarRadius, absWallW, topDBRpairs, botDBRpairs, unEtchedBotDBRpairs, minUnEtchedBotDBRpairs, maxUnEtchedBotDBRpairs, minApertureSize, apertureThickness, apertureSpacing, apertureEnable, permApertureOx, permApertureNonOx, permSubstrate, permBotDBR1, permBotDBR2, permCavity, permTopDBR1, permTopDBR2, permCap, permAbsWall, simFEMdeg, simPrecision, simFieldRes, modeNumber, jcmMaxThreads, startWL, stopWL, targetWL):
    # Set CSV field names for saving data points
    fieldnames = ['Radius [nm]','DBR detune', 'Top DBR pairs', 'unetched DBR pairs', 'Oxide width [nm]', 'nPhi', 'Polarisation', 'Wavelength [nm]', 'Q-Factor', 'Mode Volume [nm^3]', 'Purcell Factor', 'Beta Factor', 'Distance wall to mode max [nm]', 'energyIntStructure', 'energyIntAir']

    # Create folders needed to store the results, skip if they already exist (scheme: $(diameter)_nm/n_Phi$(nPhi)/*files*)
    if not os.path.exists('res_dbretchscan'):
        os.makedirs('res_dbretchscan')
        print('Creating directory \"res_dbretchscan\"')
    else:
        print('Directory \"res_dbretchscan\" already exists. Skipping creation.')
    if not os.path.exists('res_dbretchscan/fields'):
        os.makedirs('res_dbretchscan/fields')
        print('Creating directory \"res_dbretchscan/fields\"')
    else:
        print('Directory \"res_dbretchscan/fields\" already exists. Skipping creation.')

    resPath = 'res_dbretchscan/dbretchscan.csv'
    # Open csv file in write mode to store results, write header
    with open(resPath, mode='w', newline='') as csvfile:
        f = csv.DictWriter(csvfile, fieldnames=fieldnames)
        f.writeheader()
        # Iterate from min to max aperture size, start JCM solver + plot results for each
        for iterUnEtchedBotDBRpairs in range(minUnEtchedBotDBRpairs, maxUnEtchedBotDBRpairs+1, 1):
            # Variable to store number of allowed unsuccessful sim runs before cancelling the run in a row
            seqNoMode = 0
            # Temporary variable used to break out of the current nPhi loop and begin the next radius
            breakVar = False
            # Call the external function to estimate nPhi for the combination of WL and radius (cavity permeability needed for WL in medium)
            nPhi = pyWGM_modeestimation.estnPhi(minPillarRadius, stopWL, permCavity)
            # Set the initial guess to stopWL (sim starts with longest and ends with shortest WL) or call estimation function for initial guess.
            # Estimation function is based on polynomial regression, not always reasonably accurate!
            # iterGuess = 2 * speed_of_light * pi / (stopWL*1e-9)
            iterGuess = pyWGM_modeestimation.estGuessY(minPillarRadius, nPhi, permCavity)

            # Automatically derived parameters
            # substrateRadius is the radius of the substrate/air/simdomain, 1.5 times pillar radius proved to work well. Keep at at least 1000 nm for small structures.
            if minPillarRadius > 2000:
                iterSubstrateRadius = minPillarRadius + minPillarRadius/2
            else:
                iterSubstrateRadius = minPillarRadius + 1000
            # airHeight sets the height of the entire simulation area (pillar height + 300 nm)
            airHeight =  300 + substrateHeight + botDBR1Height + botDBRpairs * (botDBR1Height + botDBR2Height) + cavityHeight + topDBR1Height + topDBRpairs * (topDBR1Height + topDBR2Height) + capHeight
            # absWallH sets the height of the absorbing layer placed on the outside of the pillar - must be same height as etched pillar
            absWallH = capHeight + topDBR1Height + topDBRpairs * (topDBR1Height + topDBR2Height) + cavityHeight + (botDBRpairs - iterUnEtchedBotDBRpairs) * (botDBR1Height + botDBR2Height)
            #
            # apertureRadius determines the radius of the oxide aperture used for site-controlled QD growth. Disable aperture for size = 0.
            apertureRadius = minPillarRadius - minApertureSize
            # aperturePosition places the aperture one apertureSpacing below the cavity
            aperturePosition = - cavityHeight/2 - apertureSpacing
            # Run the  loop until the WL of the last found mode is <900 nm
            for iternPhi in range(nPhi, nPhi+20):
                # Start with the number of modes to look for from main.py. 2 usually is enough to find both WGMs for a given nPhi.
                # For very small pillars This may not be the case, as the simulation tends to lock onto higher order modes, or modes inside the DBR.
                # In this case, the simulation will be re-run with a higher number (6) of modes to look for. This will take longer, but should find 
                # both modes.
                iterModeNumber = modeNumber
                while True:
                    # Initialize variable to 0, so it reflects the number of WGMs found in one pass
                    modesFound = 0

                    # Run the JCM solver, evaluate results. Parameters passed to the solver starting with "iter" are varied between sim runs, others are fixed and taken from structure setup
                    print('Starting JCM solver for bottom DBR etch scan - unetched pairs = {} nm, nPhi = {}'.format(iterUnEtchedBotDBRpairs, iternPhi))
                    try:
                        JCMresults = pyWGM_runjcm.runJCM(jcmMaxThreads, simFEMdeg, simPrecision, simFieldRes, iterModeNumber, iterGuess, iternPhi, substrateHeight, botDBR1Height, botDBR2Height,   cavityHeight, topDBR1Height, topDBR2Height, capHeight, airHeight, absWallH, absWallW, minPillarRadius, iterSubstrateRadius, topDBRpairs, botDBRpairs,     iterUnEtchedBotDBRpairs, permSubstrate, permBotDBR1, permBotDBR2, permCavity, permTopDBR1, permTopDBR2, permCap, permAbsWall, apertureEnable, permApertureOx,   permApertureNonOx, apertureRadius, aperturePosition, apertureThickness)
                        # JCMresults: [0: wavelength [nm], 1: qFactor, 2: densityDiff, 3: vertical cross-section field strength, 4: horizontal cross-section field strength, 5: X coords in um, 6: Y    coords in um, 7: Z coords in um, 8: Re(eigenmode), 9: field-energy integrated over structure, 10: field energy integrated over air]

                        # Check if there is less than 2 modes with more field strength inside than outside of the cavity.
                        # If yes, increase the number of modes to scan for to 6, and restart loop.
                        # If it already was at 6, break loop and go to next nPhi/diameter as it is unlikely looking for more modes will find any.
                        # This break happens after the (possibly) found mode was saved.
                        for i in range(0, iterModeNumber):
                            if JCMresults[2][i] > 0 :
                                modesFound += 1
                        if modesFound < 2 and iterModeNumber < 6:
                            iterModeNumber = 6
                            print('{} whispering gallery mode(s) found. Re-scanning for {} modes instead of {}.'.format(modesFound, iterModeNumber, modeNumber))
                            continue

                        for i in range(0, iterModeNumber):
                            # Only evaluate modes where more field strength is inside the cavity than outside
                            if JCMresults[2][i] > 0:
                                # Evaluate the field
                                fields = pyWGM_runjcm.evalJCMfields(JCMresults[3][i], JCMresults[4][i], JCMresults[5], JCMresults[7], simFieldRes, JCMresults[0][i]*1e-9, JCMresults[1][i], i,  permCavity)
                                # fields: [0: field direction (radial or Y), 1:intensityXY, 2:intensityXZ, 3: x pos of field strength max [nm], 4: mode volume [nm^3], 5:purcell factor, 6: beta    factor

                                # Set plot title
                                plotTitle = 'R = %d nm, Unetched pairs = %d nm, Nphi = %d, WL = %.2f nm, Q = %.2E' %(minPillarRadius, iterUnEtchedBotDBRpairs, iternPhi, JCMresults[0][i], JCMresults   [1][i])
                                plotPath = 'res_dbretchscan/fields/{}nm_botdbrunetched_{}_{}nPhi_{}_{:.2f}nm.png'.format(minPillarRadius, iterUnEtchedBotDBRpairs, iternPhi, fields[0], JCMresults[0]   [i])
                                pyWGM_plot.plotfields(plotTitle, plotPath, JCMresults[5], JCMresults[6], JCMresults[7], fields[0], fields[1], fields[2], botDBR1Height, botDBR2Height, cavityHeight,    topDBR1Height, topDBR2Height, capHeight, topDBRpairs, botDBRpairs, iterUnEtchedBotDBRpairs, minPillarRadius, iterSubstrateRadius, apertureRadius, aperturePosition,    apertureEnable, apertureThickness)

                                # Print saved file info to the console
                                print('Plot for R = {} nm, Unetched bottom DBR pairs = {} nm, WL = {:.2f} nm saved. Mode Volume: {:.0f} nm³, purcell Factor = {:.1f}, Distance wall to mode max = {:.   1f} nm'.format(minPillarRadius, iterUnEtchedBotDBRpairs, JCMresults[0][i], fields[4], fields[5], minPillarRadius-fields[3]))
                                # Save all computed values for WGMs to the csv
                                f.writerow({'Radius [nm]': minPillarRadius, 'Top DBR pairs': topDBRpairs, 'unetched DBR pairs': iterUnEtchedBotDBRpairs, 'Oxide width [nm]':minApertureSize, 'nPhi':    iternPhi, 'Polarisation': fields[0], 'Wavelength [nm]': JCMresults[0][i], 'Q-Factor': JCMresults[1][i], 'Mode Volume [nm^3]': fields[4], 'Purcell Factor': fields  [5], 'Beta Factor': fields[6], 'Distance wall to mode max [nm]': minPillarRadius-fields[3], 'energyIntStructure': JCMresults[9][i], 'energyIntAir': JCMresults[10]    [i] })

                                # Set seqNoMode back to 0 after WGM is found
                                seqNoMode = 0
                                # Break after finding a WG mode below 900nm
                                if JCMresults[0][i] < startWL:
                                    breakVar = True

                            # If no WGM is found, increase counter by 1
                            else:
                                seqNoMode += 1 

                        # Break after either finding 2 modes, or not finding 2 even in the second pass.
                        if modesFound == 2 or iterModeNumber == 6:
                            break

                    except RuntimeError as err:
                        print(f"Fields evaluation encountered an error - Unexpected {err=}, {type(err)=}")
                    except OSError as err:
                        print(f"JCM solver encountered an error - Unexpected {err=}, {type(err)=}")
                    except BaseException as err:
                        print(f"Unknown error - Unexpected {err=}, {type(err)=}")

                # Use the current WGM circular frequency as a starting value the next time JCM is called
                # or call estimating function
                # iterGuess = JCMresults[8][i]
                iterGuess = pyWGM_modeestimation.estGuessY(minPillarRadius, iternPhi+1, permCavity)
                if breakVar:
                    break
                if seqNoMode > 12:
                    break

    # Set plot labels/paths and pass to the function for plotting result overview graphs
    plotPath= 'res_dbretchscan'
    resPath = '{}/dbretchscan.csv'.format(plotPath)
    pyWGM_plot.plotresults(resPath, plotPath, 'dbretch', 'singlePoint', targetWL)
