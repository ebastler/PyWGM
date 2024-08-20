import pyWGM_plot

# Options: 'dia', 'aperture', 'topdbr', 'dbretch', 'refindexcavity'
simType = 'dia'
# Choose the FSR determining algorithm
# average: Take the average FSR of all modes for each datapoint
# singlePoint: Take the distance of the first mode found above and below the given wavelength
# FSRtype = 'singlePoint'
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
    
    
pyWGM_plot.plotresults(resPath, plotPath, simType, FSRtype, singlepointWavelength)