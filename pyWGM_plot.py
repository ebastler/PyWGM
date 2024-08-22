import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import h,speed_of_light
import csv 
from itertools import zip_longest
import statistics
# Marker size
plotMarkerSize = 5
# How many sigma to plot
sigma = 2
# Size of the plot
plotsize_f = (8, 5)
plotsize_s = (5, 4)

def plotresults(resPath, plotPath, simType, FSRtype, singlePointWavelength):
    results = np.genfromtxt(resPath, dtype=None, delimiter=',', names=True, encoding=None)
    WL_R = []
    WL_Y = []
    plotX_R = []
    plotX_Y = []
    qFactor_Y = []
    qFactor_R = []
    modeVolume_Y=[]
    modeVolume_R=[]
    purcellF_Y=[]
    purcellF_R=[]
    betaF_Y=[]
    betaF_R=[]
    distToWall_Y=[]
    distToWall_R=[]
    FSR_Y_mode1=[]
    FSR_Y_mode2=[]
    FSR_Y_nm=[]
    FSR_Y_meV=[]
    FSR_Y_d=[]
    FSR_R_mode1=[]
    FSR_R_mode2=[]
    FSR_R_nm=[]
    FSR_R_meV=[]
    FSR_R_d=[]
    energyRatio_Y=[]
    energyRatio_R=[]
    plotX_Y_mean=[]
    plotX_R_mean=[]
    betaF_Y_mean=[]
    betaF_Y_stdev=[]
    betaF_R_mean=[]
    betaF_R_stdev=[]
    distToWall_Y_mean=[]
    distToWall_Y_stdev=[]
    distToWall_R_mean=[]
    distToWall_R_stdev=[]
    modeVolume_Y_mean=[]
    modeVolume_Y_stdev=[]
    modeVolume_R_mean=[]
    modeVolume_R_stdev=[]
    purcellF_Y_mean=[]
    purcellF_Y_stdev=[]
    purcellF_R_mean=[]
    purcellF_R_stdev=[]
    qFactor_Y_mean=[]
    qFactor_Y_stdev=[]
    qFactor_R_mean=[]
    qFactor_R_stdev=[]
    energyRatio_Y_mean=[]
    energyRatio_Y_stdev=[]
    energyRatio_R_mean=[]
    energyRatio_R_stdev=[]

    if simType == 'dia':
        plotX = (results['Radius_nm']*2)/1000
        simParameter = 'diameter'
        plotXlab = r'Pillar diameter in $\mu$m'
    if simType == 'topdbr':
        plotX = results['Top_DBR_pairs']
        simParameter = 'top DBR pairs'
        plotXlab = "Top DBR mirror pairs"
    if simType == 'refindexcavity':
        plotX = results['Refindex_Cavity']
        simParameter = 'refractive index of cavity'
        plotXlab = "Cavity refractive index"
    if simType == 'aperture':
        plotX = results['Oxide_width_nm']
        simParameter = 'oxide width'
        plotXlab = "Oxide aperture width in nm"
    if simType == 'dbretch':
        plotX = results['unetched_DBR_pairs']
        simParameter = 'unetched DBR pairs'
        plotXlab = "Unetched bottom DBR mirror pairs"


    for i in range(0, len(results)):
        # Append values to plotting arrays
        if results['Polarisation'][i] == 'Y':
            WL_Y.append(results['Wavelength_nm'][i])
            plotX_Y.append(plotX[i])
            qFactor_Y.append(results['QFactor'][i])
            modeVolume_Y.append(results['Mode_Volume_nm3'][i]*1e-9)
            purcellF_Y.append(results['Purcell_Factor'][i])
            betaF_Y.append(results['Beta_Factor'][i])
            distToWall_Y.append(results['Distance_wall_to_mode_max_nm'][i])
            energyRatio_Y.append(100*results['energyIntStructure'][i] / (results['energyIntStructure'][i] + results['energyIntAir'][i]))
            #
        else:
            WL_R.append(results['Wavelength_nm'][i])
            plotX_R.append(plotX[i])
            qFactor_R.append(results['QFactor'][i])
            modeVolume_R.append(results['Mode_Volume_nm3'][i]*1e-9)
            purcellF_R.append(results['Purcell_Factor'][i])
            betaF_R.append(results['Beta_Factor'][i])
            distToWall_R.append(results['Distance_wall_to_mode_max_nm'][i])
            energyRatio_R.append(100*results['energyIntStructure'][i] / (results['energyIntStructure'][i] + results['energyIntAir'][i]))



    if FSRtype == 'singlePoint':
    # Determine FSR for singlePoint method
    # Keep scanning data as long as WL > singlePoint and store the last checked i value to a temp variable
    # Once WL is smaller than singlePoint, check if tempvar != 0 (last one was already smaller than singlePoint -> 0) and x value is the same.
    # Store difference between temp variable (last WL > singlePoint) and current variable (first WL < singlePoint) to file, then set temp variable to 0
    # Will remain 0 until WL > singlePoint (new x value).
        tempvarPrevY = -1
        tempvarPrevR = -1
        for i in range(0, len(results)-1):
            if results['Polarisation'][i] == 'Y':
                if results['Wavelength_nm'][i] > singlePointWavelength:
                    tempvarPrevY = i
                else:
                    if tempvarPrevY >= 0:
                        if plotX[i] == plotX[tempvarPrevY]:
                            FSR_Y_mode1.append(results['Wavelength_nm'][i])
                            FSR_Y_mode2.append(results['Wavelength_nm'][tempvarPrevY])
                            FSR_Y_nm.append(abs(results['Wavelength_nm'][i] - results['Wavelength_nm'][tempvarPrevY]))
                            FSR_Y_meV.append(abs(h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][i] - h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][tempvarPrevY]))
                            FSR_Y_d.append(plotX[i])
                            tempvarPrevY = -1

            if results['Polarisation'][i] == 'radial':
                if results['Wavelength_nm'][i] > singlePointWavelength:
                    tempvarPrevR = i
                else:
                    if tempvarPrevR >= 0:
                        if plotX[i] == plotX[tempvarPrevR]:
                            FSR_R_mode1.append(results['Wavelength_nm'][i])
                            FSR_R_mode2.append(results['Wavelength_nm'][tempvarPrevR])
                            FSR_R_nm.append(abs(results['Wavelength_nm'][i] - results['Wavelength_nm'][tempvarPrevR]))
                            FSR_R_meV.append(abs(h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][i] - h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][tempvarPrevR]))
                            FSR_R_d.append(plotX[i])
                            tempvarPrevR = -1
        
    if FSRtype == 'fullData':
    # Determine FSR for fullData method
    # Compute FSR between each 2 neighboring modes of the same diameter
    # Once WL is smaller than singlePoint, check if tempvar != 0 (last one was already smaller than singlePoint -> 0) and x value is the same.
    # Store difference between temp variable (last WL > singlePoint) and current variable (first WL < singlePoint) to file, then set temp variable to 0
    # Will remain 0 until WL > singlePoint (new x value).
        tempvarPrevY = -1
        tempvarPrevR = -1
        for i in range(0, len(results)-1):
            print()
            if results['Polarisation'][i] == 'Y':
                if tempvarPrevY >= 0:
                    if plotX[i] == plotX[tempvarPrevY]:
                        FSR_Y_mode1.append(results['Wavelength_nm'][i])
                        FSR_Y_mode2.append(results['Wavelength_nm'][tempvarPrevY])
                        FSR_Y_nm.append(abs(results['Wavelength_nm'][i] - results['Wavelength_nm'][tempvarPrevY]))
                        FSR_Y_meV.append(abs(h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][i] - h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][tempvarPrevY]))
                        FSR_Y_d.append(plotX[i])
                tempvarPrevY = i

            if results['Polarisation'][i] == 'radial':
                if tempvarPrevR >= 0:
                    if plotX[i] == plotX[tempvarPrevR]:
                        FSR_R_mode1.append(results['Wavelength_nm'][i])
                        FSR_R_mode2.append(results['Wavelength_nm'][tempvarPrevR])
                        FSR_R_nm.append(abs(results['Wavelength_nm'][i] - results['Wavelength_nm'][tempvarPrevR]))
                        FSR_R_meV.append(abs(h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][i] - h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][tempvarPrevR]))
                        FSR_R_d.append(plotX[i])
                tempvarPrevR = i

                
    if FSRtype == 'average':
    # Determine FSR for average method
    # If line i is a Y polarized mode, check for next Y polarized line j, append the abs(distance) between them to tempvar.
    # Break if the radius changes or a mode is found.
        tempvar_nm = []
        tempvar_meV = []
        for i in range(0, len(results)-1):
            if results['Polarisation'][i] == 'Y':
                for j in range(i+1, len(results)-1):
                    if plotX[i] == plotX[j]:
                        if results['Polarisation'][j] == 'Y':
                            tempvar_nm.append(abs(results['Wavelength_nm'][i] - results['Wavelength_nm'][j]))
                            tempvar_meV.append(abs(h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][i] - h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][j]))
                            break
            if plotX[i] != plotX[i+1]:
                FSR_Y_nm.append(np.mean(tempvar_nm))
                FSR_Y_meV.append(np.mean(tempvar_meV))
                FSR_Y_d.append(plotX[i])
                tempvar_nm = []
                tempvar_meV = []
        FSR_Y_nm.append(np.mean(tempvar_nm))
        FSR_Y_meV.append(np.mean(tempvar_meV))
        FSR_Y_d.append(plotX[i])

        tempvar_nm = []
        tempvar_meV = []
        for i in range(0, len(results)-1):
            if results['Polarisation'][i] == 'radial':
                for j in range(i+1, len(results)-1):
                    if plotX[i] == plotX[j]:
                        if results['Polarisation'][j] == 'radial':
                            tempvar_nm.append(abs(results['Wavelength_nm'][i] - results['Wavelength_nm'][j]))
                            tempvar_meV.append(abs(h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][i] - h*speed_of_light/1.602176634E-19/1E-12/results['Wavelength_nm'][j]))
                            break
            if plotX[i] != plotX[i+1]:
                FSR_R_nm.append(np.mean(tempvar_nm))
                FSR_R_meV.append(np.mean(tempvar_meV))
                FSR_R_d.append(plotX[i])
                tempvar_nm = []
                tempvar_meV = []
        FSR_R_nm.append(np.mean(tempvar_nm))
        FSR_R_meV.append(np.mean(tempvar_meV))
        FSR_R_d.append(plotX[i])

    # Write FSR data to save files
    data = [FSR_Y_mode1, FSR_Y_mode2, FSR_Y_d, FSR_Y_meV, FSR_Y_nm]
    export_data = zip_longest(*data, fillvalue = '')
    with open('{}/fsr_y.csv'.format(plotPath), 'w', encoding="ISO-8859-1", newline='') as file:
        write = csv.writer(file)
        write.writerow(("Mode 1", "Mode 2", "FSR_Y_d", "FSR_Y_meV", "FSR_Y_nm"))
        write.writerows(export_data)

    data = [FSR_R_mode1, FSR_R_mode2, FSR_R_d, FSR_R_meV, FSR_R_nm]
    export_data = zip_longest(*data, fillvalue = '')
    with open('{}/fsr_r.csv'.format(plotPath), 'w', encoding="ISO-8859-1", newline='') as file:
        write = csv.writer(file)
        write.writerow(("Mode 1", "Mode 2", "FSR_Y_d", "FSR_Y_meV", "FSR_Y_nm"))
        write.writerows(export_data)

    # Calculate mean of beta-factor etc for all simulated modes per wavelength (i.e. average/stdev for all vertically polarized modes in 1000nm pillar, etc)
    tempArrayY = []
    tempArrayR = []
    # Loop over all values for vertical polarization (Y)
    for i in range(0, len(plotX_Y)):
        # If array 0 (first value) or the previous diameter was the same as the current diameter, append i to temp array
        if len(tempArrayY) == 0 or plotX_Y[tempArrayY[-1]] == plotX_Y[i]:
            tempArrayY.append(i)
        # At diameter change, use the values for all previously determined indizes for mean/stdev
        else:
            betaF_Y_mean.append(statistics.mean([betaF_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
            betaF_Y_stdev.append(statistics.stdev([betaF_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
            distToWall_Y_mean.append(statistics.mean([distToWall_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
            distToWall_Y_stdev.append(statistics.stdev([distToWall_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
            modeVolume_Y_mean.append(statistics.mean([modeVolume_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
            modeVolume_Y_stdev.append(statistics.stdev([modeVolume_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
            purcellF_Y_mean.append(statistics.mean([purcellF_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
            purcellF_Y_stdev.append(statistics.stdev([purcellF_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
            qFactor_Y_mean.append(statistics.mean([qFactor_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
            qFactor_Y_stdev.append(statistics.stdev([qFactor_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
            energyRatio_Y_mean.append(statistics.mean([energyRatio_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
            energyRatio_Y_stdev.append(statistics.stdev([energyRatio_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
            plotX_Y_mean.append(plotX_Y[i-1])
            tempArrayY = []
            tempArrayY.append(i)
    # Append last values to arrays because loop above only appends at diameter change
    betaF_Y_mean.append(statistics.mean([betaF_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
    betaF_Y_stdev.append(statistics.stdev([betaF_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
    distToWall_Y_mean.append(statistics.mean([distToWall_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
    distToWall_Y_stdev.append(statistics.stdev([distToWall_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
    modeVolume_Y_mean.append(statistics.mean([modeVolume_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
    modeVolume_Y_stdev.append(statistics.stdev([modeVolume_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
    purcellF_Y_mean.append(statistics.mean([purcellF_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
    purcellF_Y_stdev.append(statistics.stdev([purcellF_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
    qFactor_Y_mean.append(statistics.mean([qFactor_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
    qFactor_Y_stdev.append(statistics.stdev([qFactor_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
    energyRatio_Y_mean.append(statistics.mean([energyRatio_Y[j] for j in tempArrayY]) if len(tempArrayY) >0 else np.nan)
    energyRatio_Y_stdev.append(statistics.stdev([energyRatio_Y[j] for j in tempArrayY]) if len(tempArrayY) >1 else 0)
    plotX_Y_mean.append(plotX_Y[i-1])

    # Loop over all values for vertical polarization (Y)
    for i in range(0, len(plotX_R)):
        # If array 0 (first value) or the previous diameter was the same as the current diameter, append i to temp array
        if len(tempArrayR) == 0 or plotX_R[tempArrayR[-1]] == plotX_R[i]:
            tempArrayR.append(i)
        # At diameter change, use the values for all previously determined indizes for mean/stdev
        else:
            betaF_R_mean.append(statistics.mean([betaF_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
            betaF_R_stdev.append(statistics.stdev([betaF_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
            distToWall_R_mean.append(statistics.mean([distToWall_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
            distToWall_R_stdev.append(statistics.stdev([distToWall_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
            modeVolume_R_mean.append(statistics.mean([modeVolume_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
            modeVolume_R_stdev.append(statistics.stdev([modeVolume_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
            purcellF_R_mean.append(statistics.mean([purcellF_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
            purcellF_R_stdev.append(statistics.stdev([purcellF_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
            qFactor_R_mean.append(statistics.mean([qFactor_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
            qFactor_R_stdev.append(statistics.stdev([qFactor_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
            energyRatio_R_mean.append(statistics.mean([energyRatio_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
            energyRatio_R_stdev.append(statistics.stdev([energyRatio_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
            plotX_R_mean.append(plotX_R[i-1])
            tempArrayR = []
            tempArrayR.append(i)
    # Append last values to arrays because loop above only appends at diameter change
    betaF_R_mean.append(statistics.mean([betaF_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
    betaF_R_stdev.append(statistics.stdev([betaF_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
    distToWall_R_mean.append(statistics.mean([distToWall_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
    distToWall_R_stdev.append(statistics.stdev([distToWall_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
    modeVolume_R_mean.append(statistics.mean([modeVolume_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
    modeVolume_R_stdev.append(statistics.stdev([modeVolume_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
    purcellF_R_mean.append(statistics.mean([purcellF_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
    purcellF_R_stdev.append(statistics.stdev([purcellF_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
    qFactor_R_mean.append(statistics.mean([qFactor_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
    qFactor_R_stdev.append(statistics.stdev([qFactor_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
    energyRatio_R_mean.append(statistics.mean([energyRatio_R[j] for j in tempArrayR]) if len(tempArrayR) >0 else np.nan)
    energyRatio_R_stdev.append(statistics.stdev([energyRatio_R[j] for j in tempArrayR]) if len(tempArrayR) >1 else 0)
    plotX_R_mean.append(plotX_R[i-1])

    # Write mean/stdev data to save files
    data = [plotX_Y_mean, betaF_Y_mean, betaF_Y_stdev, distToWall_Y_mean, distToWall_Y_stdev, modeVolume_Y_mean, modeVolume_Y_stdev, purcellF_Y_mean, purcellF_Y_stdev, qFactor_Y_mean, qFactor_Y_stdev, energyRatio_Y_mean, energyRatio_Y_stdev]
    export_data = zip_longest(*data, fillvalue = '')
    with open('{}/meanvalues_y.csv'.format(plotPath), 'w', encoding="ISO-8859-1", newline='') as file:
        write = csv.writer(file)
        write.writerow(("Diameter_nm", "BetaFactor_mean", "BetaFactor_stdev", "DistToWall_mean", "DistToWall_stdev", "ModeVolume_mean", "ModeVolume_stdev", "PurcellFactor_mean", "PurcellFactor_stdev", "qFactor_mean", "qFactor_stdev", "EnergyRatio_mean", "EnergyRatio_stdev"))
        write.writerows(export_data)

    data = [plotX_R_mean, betaF_R_mean, betaF_R_stdev, distToWall_R_mean, distToWall_R_stdev, modeVolume_R_mean, modeVolume_R_stdev, purcellF_R_mean, purcellF_R_stdev, qFactor_R_mean, qFactor_R_stdev, energyRatio_R_mean, energyRatio_R_stdev]
    export_data = zip_longest(*data, fillvalue = '')
    with open('{}/meanvalues_r.csv'.format(plotPath), 'w', encoding="ISO-8859-1", newline='') as file:
        write = csv.writer(file)
        write.writerow(("Radius_nm", "BetaFactor_mean", "BetaFactor_stdev", "DistToWall_mean", "DistToWall_stdev", "ModeVolume_mean", "ModeVolume_stdev", "PurcellFactor_mean", "PurcellFactor_stdev", "qFactor_mean", "qFactor_stdev", "EnergyRatio_mean", "EnergyRatio_stdev"))
        write.writerows(export_data)
            
    plotTitle = ('FSR and Q-factor for dia = %d nm to dia = %d nm'  %(min(results['Radius_nm'])*2, max(results['Radius_nm'])*2))

    plotgraphsingle(simParameter, plotTitle, plotPath, plotXlab, plotX_R, WL_R, qFactor_R, modeVolume_R, purcellF_R, betaF_R, distToWall_R, plotX_Y, WL_Y, qFactor_Y, modeVolume_Y, purcellF_Y, betaF_Y, distToWall_Y, FSRtype, singlePointWavelength, FSR_Y_d, FSR_Y_nm, FSR_Y_meV, FSR_R_d, FSR_R_nm, FSR_R_meV, energyRatio_Y, energyRatio_R, plotX_Y_mean, betaF_Y_mean, [n*sigma for n in betaF_Y_stdev] , distToWall_Y_mean, [n*sigma for n in distToWall_Y_stdev], modeVolume_Y_mean, [n*sigma for n in modeVolume_Y_stdev], purcellF_Y_mean, [n*sigma for n in purcellF_Y_stdev], qFactor_Y_mean, [n*sigma for n in qFactor_Y_stdev], energyRatio_Y_mean, [n*sigma for n in energyRatio_Y_stdev], plotX_R_mean, betaF_R_mean, [n*sigma for n in betaF_R_stdev], distToWall_R_mean, [n*sigma for n in distToWall_R_stdev], modeVolume_R_mean, [n*sigma for n in modeVolume_R_stdev], purcellF_R_mean, [n*sigma for n in purcellF_R_stdev], qFactor_R_mean, [n*sigma for n in qFactor_R_stdev], energyRatio_R_mean, [n*sigma for n in energyRatio_R_stdev], plotsize_f, plotsize_s)

def plotfields(plotTitle, plotPath, posX, posY, posZ, intensityXZ_dir, intensityXY, intensityXZ, botDBR1Height, botDBR2Height, cavityHeight, topDBR1Height, topDBR2Height, capHeight, topDBRpairs, botDBRpairs, unEtchedBotDBRpairs, pillarRadius, substrateRadius, apertureRadius, aperturePosition, apertureEnable, apertureThickness):
    plt.rc('font', family='serif', size=14)
    matplotlib.rc('text', usetex=True)
    # Plot the results as a heatmap, save as png
    fig, axes = plt.subplots(1,2, figsize=(10, 4))
    #
    # Plot title
    plt.suptitle(plotTitle)
    #
    # Main plot, vertical cross-section
    im = axes[0].pcolormesh(posX, posY, intensityXY, vmin=-1., vmax=1., shading='gouraud', cmap='seismic')
    # Add pillar outlines
    pillarTop = (cavityHeight/2 + topDBR1Height + topDBRpairs*(topDBR1Height + topDBR2Height) + capHeight)/1000  
    pillarBase = (-cavityHeight/2 - (botDBRpairs - unEtchedBotDBRpairs) * (botDBR1Height + botDBR2Height))/1000
    cavityTop = cavityHeight/2000
    pillarR = pillarRadius/1000
    substrateR = substrateRadius/1000
    substrateH = (-cavityHeight/2 - botDBR1Height - botDBRpairs * (botDBR1Height + botDBR2Height))/1000
    axes[0].plot([-pillarR, pillarR], [cavityTop, cavityTop], 'k-', lw=1)
    axes[0].plot([-pillarR, pillarR], [-cavityTop, -cavityTop], 'k-', lw=1)
    axes[0].plot([-pillarR, pillarR], [pillarTop, pillarTop], 'k-', lw=1)
    axes[0].plot([-pillarR, -substrateR], [pillarBase, pillarBase], 'k-', lw=1)
    axes[0].plot([pillarR, substrateR], [pillarBase, pillarBase], 'k-', lw=1)
    axes[0].plot([-substrateR, substrateR], [substrateH, substrateH], 'k-', lw=1)
    axes[0].plot([-pillarR, -pillarR], [pillarBase, pillarTop], 'k-', lw=1)
    axes[0].plot([pillarR, pillarR], [pillarBase, pillarTop], 'k-', lw=1)
    if apertureEnable:
        apertureR = apertureRadius/1000
        apertureH = aperturePosition/1000
        apertureT = apertureThickness/1000
        axes[0].plot([-pillarR, -apertureR, -apertureR, -pillarR], [apertureH, apertureH, apertureH-apertureT, apertureH-apertureT], 'g-', lw=1)
        axes[0].plot([pillarR, apertureR, apertureR, pillarR], [apertureH, apertureH, apertureH-apertureT, apertureH-apertureT], 'g-', lw=1)
    # Force all axes to have the same scale
    axes[0].set_aspect('equal', 'box')
    # Add x and y axis label
    axes[0].set(xlabel = r'Pillar width (x) in $\mu$m', ylabel = r'Pillar height (y) in $\mu$m')
    # Add colorbar (Z axis legend)
    fig.colorbar(im, ax=axes[0], label='Field strength in %s direction (a.u.)' %(intensityXZ_dir))
    #
    # Main plot, horizontal cross-section
    im = axes[1].pcolormesh(posX, posZ, intensityXZ, vmin=-1., vmax=1., shading='gouraud', cmap='seismic')
    # Add circle with pillar diameter
    circle1 = plt.Circle((0,0), pillarR, ec='black', fill=False)
    axes[1].add_patch(circle1)
    # Force all axes to have the same scale
    axes[1].set_aspect('equal', 'box')
    # Add x and y axis label
    axes[1].set(xlabel = r'Pillar width (x) in $\mu$m', ylabel = r'Pillar width (z) in $\mu$m')
    # Add colorbar (Z axis legend)
    fig.colorbar(im, ax=axes[1], label='Field strength in %s direction (a.u.)' %(intensityXZ_dir))
    #
    # Save plots to file
    plt.tight_layout()
    plt.savefig(plotPath, dpi=200)
    plt.close('all')

def plotgraphsingle(simParameter, plotTitle, plotPath, plotXlab, plotXZ_X, plotXZ_WL, plotXZ_q, plotXZ_mv, plotXZ_pf, plotXZ_bf, plotXZ_dist, plotY_X, plotY_WL, plotY_q, plotY_mv, plotY_pf, plotY_bf, plotY_dist, FSRtype, singlePointWavelength, FSR_Y_d, FSR_Y_nm, FSR_Y_meV, FSR_R_d, FSR_R_nm, FSR_R_meV, energyRatio_Y, energyRatio_R, plotX_Y_mean, betaF_Y_mean, betaF_Y_stdev, distToWall_Y_mean, distToWall_Y_stdev, modeVolume_Y_mean, modeVolume_Y_stdev, purcellF_Y_mean, purcellF_Y_stdev, qFactor_Y_mean, qFactor_Y_stdev, energyRatio_Y_mean, energyRatio_Y_stdev, plotX_R_mean, betaF_R_mean, betaF_R_stdev, distToWall_R_mean, distToWall_R_stdev, modeVolume_R_mean, modeVolume_R_stdev, purcellF_R_mean, purcellF_R_stdev, qFactor_R_mean, qFactor_R_stdev, energyRatio_R_mean, energyRatio_R_stdev, plotsize_f, plotsize_s):
    #
    plt.rc('font', family='serif', size=14)
    matplotlib.rc('text', usetex=True)
    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Wavelength over {}'.format(simParameter))
    ax1.scatter(plotY_X, plotY_WL, s=plotMarkerSize, label = 'Vertically polarized mode')
    ax1.scatter(plotXZ_X, plotXZ_WL, s=plotMarkerSize, label = 'Radially polarized mode')
    ax1.set(xlabel = plotXlab, ylabel = "Wavelength in nm")
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_WL.pdf'.format(plotPath))
    plt.close()
    #
    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #ax2 = ax1.twinx()
    #plt.suptitle('FSR over {}'.format(simParameter))
    ax1.scatter(FSR_Y_d, FSR_Y_nm, s=plotMarkerSize, label = 'Vertically polarized mode')
    ax1.scatter(FSR_R_d, FSR_R_nm, s=plotMarkerSize, label = 'Radially polarized mode')
    #ax2.scatter(FSR_Y_d, FSR_Y_nm, s=0.5)
    #ax2.scatter(FSR_R_d, FSR_R_nm, s=0.5)
    if FSRtype == 'singlePoint':
        ax1.set(xlabel = plotXlab, ylabel = "FSR at {} nm in nm".format(singlePointWavelength))
    if FSRtype == 'average':
        ax1.set(xlabel = plotXlab, ylabel = "Mean FSR in nm")
    if FSRtype == 'fullData':
        ax1.set(xlabel = plotXlab, ylabel = "FSR in nm")
    
    #ax2.set(ylabel = "Mean FSR in nm")
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_FSR.pdf'.format(plotPath))
    plt.close()
    #
    # Plot means
    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Mean beta-factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, betaF_Y_mean, yerr=betaF_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, betaF_R_mean, yerr=betaF_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = r'$\beta$-factor')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_mean_Beta.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Distance from mode max to wall over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, distToWall_Y_mean, yerr=distToWall_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, distToWall_R_mean, yerr=distToWall_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = "Distance mode max to wall in nm")
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_mean_Dist.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Mode volume over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, modeVolume_Y_mean, yerr=modeVolume_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, modeVolume_R_mean, yerr=modeVolume_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = r'Mode Volume in $\mu$m³')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_mean_MV.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Purcell factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, purcellF_Y_mean, yerr=purcellF_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, purcellF_R_mean, yerr=purcellF_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = "Purcell factor")
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_mean_Purcell.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Purcell factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, purcellF_Y_mean, yerr=purcellF_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, purcellF_R_mean, yerr=purcellF_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = "Purcell factor")
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_mean_Purcell.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Q-factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, qFactor_Y_mean, yerr=qFactor_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, qFactor_R_mean, yerr=qFactor_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = "Q-factor")
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_mean_Q.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_f)
    #plt.suptitle('Mode energy inside/outside structure over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, energyRatio_Y_mean, yerr=energyRatio_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, energyRatio_R_mean, yerr=energyRatio_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = r'\% of mode energy inside structure')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_f_mean_Ratio.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Wavelength over {}'.format(simParameter))
    ax1.scatter(plotY_X, plotY_WL, s=plotMarkerSize, label = 'Vertically polarized mode')
    ax1.scatter(plotXZ_X, plotXZ_WL, s=plotMarkerSize, label = 'Radially polarized mode')
    ax1.set(xlabel = plotXlab, ylabel = 'Wavelength in nm')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_WL.pdf'.format(plotPath))
    plt.close()
    #
    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #ax2 = ax1.twinx()
    #plt.suptitle('FSR over {}'.format(simParameter))
    ax1.scatter(FSR_Y_d, FSR_Y_nm, s=plotMarkerSize, label = 'Vertically polarized mode')
    ax1.scatter(FSR_R_d, FSR_R_nm, s=plotMarkerSize, label = 'Radially polarized mode')
    #ax2.scatter(FSR_Y_d, FSR_Y_nm, s=0.5)
    #ax2.scatter(FSR_R_d, FSR_R_nm, s=0.5)
    if FSRtype == 'singlePoint':
        ax1.set(xlabel = plotXlab, ylabel = "FSR at {} nm in nm".format(singlePointWavelength))
    if FSRtype == 'average':
        ax1.set(xlabel = plotXlab, ylabel = "Mean FSR in nm")
    if FSRtype == 'fullData':
        ax1.set(xlabel = plotXlab, ylabel = "FSR in nm")

    #ax2.set(ylabel = 'Mean FSR in nm')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_FSR.pdf'.format(plotPath))
    plt.close()
    #
    # Plot means
    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Mean beta-factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, betaF_Y_mean, yerr=betaF_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, betaF_R_mean, yerr=betaF_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = r'$\beta$-factor')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_mean_Beta.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Distance from mode max to wall over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, distToWall_Y_mean, yerr=distToWall_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, distToWall_R_mean, yerr=distToWall_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = 'Distance mode max to wall in nm')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_mean_Dist.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Mode volume over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, modeVolume_Y_mean, yerr=modeVolume_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, modeVolume_R_mean, yerr=modeVolume_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = r'Mode Volume in $\mu$m³')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_mean_MV.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Purcell factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, purcellF_Y_mean, yerr=purcellF_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, purcellF_R_mean, yerr=purcellF_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = 'Purcell factor')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_mean_Purcell.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Purcell factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, purcellF_Y_mean, yerr=purcellF_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, purcellF_R_mean, yerr=purcellF_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = 'Purcell factor')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_mean_Purcell.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Q-factor over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, qFactor_Y_mean, yerr=qFactor_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, qFactor_R_mean, yerr=qFactor_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = 'Q-factor')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_mean_Q.pdf'.format(plotPath))
    plt.close()

    fig, ax1 = plt.subplots(1,1, figsize=plotsize_s)
    #plt.suptitle('Mode energy inside/outside structure over {}'.format(simParameter))
    ax1.errorbar(plotX_Y_mean, energyRatio_Y_mean, yerr=energyRatio_Y_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Vertically polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.errorbar(plotX_R_mean, energyRatio_R_mean, yerr=energyRatio_R_stdev, fmt='o', ms=plotMarkerSize/2, capsize=5, label = r'Radially polarized (${}\sigma$ intervals)'.format(sigma))
    ax1.set(xlabel = plotXlab, ylabel = r'\% of mode energy inside structure')
    plt.tight_layout()
    ax1.legend()
    plt.savefig('{}/plot_s_mean_Ratio.pdf'.format(plotPath))
    plt.close()