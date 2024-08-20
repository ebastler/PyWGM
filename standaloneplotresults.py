import pyWGM_plot

# Options:  'aperture', 'dia', 'dbretch', 'topdbr', 'refindexcavity'
simType = 'dia'
# Choose the FSR determining algorithm
# average: Take the average FSR of all modes for each datapoint
# singlePoint: Take the distance of the first mode found above and below the given wavelength
# FSRtype = 'fullData'
# FSRtype = 'average'
FSRtype = 'singlePoint'
singlepointWavelength = 960

if simType == 'refindexcavity':
    plotPath= 'res_refindexscan'
    resPath = '{}/refindexscan.csv'.format(plotPath)
if simType == 'dia':
    plotPath= 'res_diascan'
    resPath = '{}/diascan.csv'.format(plotPath)
if simType == 'topdbr':
    plotPath= 'res_topdbrscan'
    resPath = '{}/topdbrscan.csv'.format(plotPath)
if simType == 'aperture':
    plotPath= 'res_aperturescan'
    resPath = '{}/aperturescan.csv'.format(plotPath)
if simType == 'dbretch':
    plotPath= 'res_dbretchscan'
    resPath = '{}/dbretchscan.csv'.format(plotPath)
    
print("Re-drawing plots for {} simulation in directory './{}'.".format(simType, plotPath))
if FSRtype == 'singlePoint':
    print("FSR type: {}. FSR wavelength: {} nm\n".format(FSRtype, singlepointWavelength))
else:
    print("FSR type: {}.\n".format(FSRtype))
    
pyWGM_plot.plotresults(resPath, plotPath, simType, FSRtype, singlepointWavelength)