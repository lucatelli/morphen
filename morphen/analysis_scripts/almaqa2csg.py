# ALMA Quality Assurance Software
# QA2 Calibration Script Generator
# Dirk Petry (ESO)
# Todd Hunter (NRAO)
# Luke Maud (ESO)
#
# $Id: almaqa2csg.py,v 2.18 2025/03/05 08:14:59 dpetry Exp $
#
"""
The ALMA QA2 Calibration Script Generator
"""

from __future__ import print_function

import os
import sys
import numpy as np
import glob
import re
import time as timeUtilities
import pprint
import analysisUtils as aU
sfsdr = aU.stuffForScienceDataReduction()

try:  # Python 3
    import inspect
    from casatasks import importasdm
    from casatasks import gencal
    from casatasks import casalog
    from casatasks import mstransform
    from casatools import table as tbtool
    from casatools import msmetadata as msmdtool
    from casatools import quanta as qatool

    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request, HTTPPasswordMgrWithDefaultRealm, HTTPBasicAuthHandler, build_opener
    
    from urllib.error import HTTPError
    from subprocess import getstatusoutput, getoutput
    import XmlObjectifier_python3 as XmlObjectifier

    def get_default_args(func):
        signature = inspect.signature(func)
        return {
            k: v.default
            for k, v in signature.parameters.items()
            if v.default is not inspect.Parameter.empty
        }

except ImportError:  # Python 2
    from taskinit import *
    from importasdm_cli import importasdm_cli as importasdm
    from gencal_cli import gencal_cli as gencal

    import exceptions

    import XmlObjectifier
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError, HTTPPasswordMgrWithDefaultRealm, HTTPBasicAuthHandler, build_opener
    from commands import getstatusoutput, getoutput


def version(short=False):
    """
    Returns the CVS revision number.
    """
    myversion = "$Id: almaqa2csg.py,v 2.18 2025/03/05 08:14:59 dpetry Exp $"
    if (short):
        myversion = myversion.split()[2]
    return myversion



def generateReducScript(msNames='', step='calib', corrAntPos=True, timeBinForFinalData=0., 
                        refant='', bpassCalId='', chanWid=1, angScale=0, run=False, lowSNR=False, 
                        projectCode='', schedblockName='', schedblockUid='', queue='', state='', 
                        upToTimeForState=2, useLocalAlmaHelper=True, tsysChanTol=1, sdQSOflux=1, 
                        runPhaseClosure=False, skipSyscalChecks=False, lazy=False, lbc=False, 
                        remcloud=False, bdfflags=True,
                        tsysPerField=False, splitMyScienceSpw=True, bpassCalTableName='', 
                        reindexMyScienceSpw=False, useCalibratorService=True,
                        calibratorServiceURL=None, allowHybrid=False,
                        combineB2BLFspws=False, combineB2BLFHFspws=False, combineB2BDGCspws=False,
                        includeRenorm=True, vetoBPBootstrap=False, legacyfixes=False):
    """
    The ALMA QA2 calibration script generator

    msNames: a string or a list of strings of UIDs (either ASDM or MS) to process 
             NOTE: rigorous regression testing is presently only done on single UIDs, not lists
              default=''
    step:    calib, fluxcal, wvr, calsurvey, SDeff, SDcalibLine, SDcalibCont, SDscience, SDampcal
              default='calib'
    corrAntPos: if True, then run correctMyAntennaPositions
              default=True
    timeBinForFinalData: a value in seconds (string, int, or float), passed to split
              default=0.
    refant:  the reference antenna to use (instead of automatic selection), must be a string
              default='', i.e. determine automatically
    bpassCalId: use the specified source for bandpass (rather than determine from the intents)
              default='', i.e. determine from the intents
    chanWid: integer, used by runCleanOnSource and searchForLines
              default=1
    angScale: deprecated

    run:     deprecated

    lowSNR:  Boolean passed to doBandpassCalibration to use whole spw for pre-bandpass phase-up
              default=False
    projectCode, schedblockName, queue, state, upToTimeForState: deprecated

    useLocalAlmaHelper: if True, run tsysspwmap inside generator, rather than in the resulting script
              default=True
    tsysChanTol: integer argument passed to tsysspwmap
              default=1
    sdQSOflux: flux density to use for quasar in single dish case (step='SDeff')
              default=1
    runPhaseClosure: deprecated

    skipSyscalChecks: if True, then don't check for negative Tsys problems
              default=False
    lazy: value of the 'lazy' parameter in importasdm. If True, reference the ASDM instead
            of copying the visibilities into the DATA column of the MS. Saves disk space.
              default=False
    lbc:    if True,  in bandpass calibration, use solint='inf,8MHz' instead of 'inf,20ch'
              default=False
    remcloud: if True, run the recipe remove_cloud prior to running wvrgcal
              default=False
    bdfflags: passed to importasdm to invoke the application of BDF flags
              default=True
    bpassCalTableName: to use instead of default name
              default='', i.e. use the bp table created for bpassCalId with the standard naming
    phaseDiff: deprecated (BWSW and B2B modes are recognized automatically)

    tsysPerField: passed to the perField parameter of tsysspwmap
              default=False
    splitMyScienceSpw: In the final split-out, only include the SPWs corresponding to intent OBSERVE_TARGET
              and BANDPASS.
              default=True
    reindexMyScienceSpw: perform reindexing in the split out after the apriori calibration.
              default=False
    useCalibratorService: if True, then, in the call to aU.getALMAFluxForMS in the setjy step, use
              aU.calibratorService(), otherwise use aU.getALMAFlux()
              default=True
    calibratorServiceURL: the URL to pass to aU.calibratorService() if useCalibratorService==True
              default: None - use the default of calibratorServiceURL in aU.getALMAFluxForMS()
    allowHybrid: if False, only the antennas of the dominant (most often occuring) antenna diameter are split out.
              If True, all antennas are split out. 
              default=False
    combineB2BLFspws: if True, and if the dataset uses band-to-band phase transfer, then combine the LF SPWs
              in gaincal for PHASE calibrator and use appropriate spwmaps. This is "CASE B", minor combine.
              default=False
    combineB2BLFHFspws: if True, and if the dataset uses band-to-band phase transfer, then combine the LF SPWs
              and HF SPWs in BANDPASS, DIFFGAIN (phaseint) and PHASE gaincal and use appropriate spwmaps. MOST USEFUL,
              required for narrow Bandwidth SPW, particualrly at HF. This is "CASE C", extension of "CASE B".
              default=False
    combineB2BDGCspws: if True, and if the dataset uses band-to-band phase transfer, then combine the LF SPWs,
              HF SPWs, and SpWs for the B2B offset in BANDPASS, DIFFGAIN, DIFFGAIN(B2B offset) and PHASE
              gaincal stages and use appropriate spwmaps. This is "CASE D", extension of "CASE C" with the added
              combine needed for the _rare_ case of weak DIFFGAIN and B2B offset with low SNR.
              default=False
    includeRenorm: if True, a step is added in the calibration script to perform renormalization
              default=True
    vetoBPBootstrap: if True, bandpass bootstrapping will not be done even if ampcal and bandpass are different
              default=False
    legacyfixes: if True, generate code to apply fixes for CSV2555 and SYSCAL table times, and fixplanets even for
              data observed on or after 1 October 2015 (ALMA Cycle 3 start). For data observed before this time,
              the fixes are applied in any case.
              default=False

    """

    print("The ALMA QA2 calibration script generator")
    print(version())
    print('using '+aU.version())
    casalog.post(version(),'INFO')
    casalog.post('using '+aU.version(),'INFO')

    mycasaversion = aU.getCasaVersion()

    if mycasaversion < '5.6.1':
        casalog.post('CASA versions < 5.6.1 are no longer supported.', 'SEVERE')

    latestValidatedVersion = '6.6.1' # 2024-10-23
    if re.search('^'+latestValidatedVersion, mycasaversion) == None:
        print('WARNING: You are currently running CASA %s rather than CASA %s.' % (mycasaversion, latestValidatedVersion))
        print('WARNING: If you observe any issue, please file an ALMA PRTSPR or helpdesk ticket.')

    if useCalibratorService:
        if calibratorServiceURL==None:
            if mycasaversion > '5.9.9':
                calibratorServiceURL = get_default_args(aU.getALMAFluxForMS)['calibratorServiceURL']
            else:
                casalog.post('Cannot determine default for calibratorServiceURL in CASA versions below 6.\nSet useCalibratorService to False or provide URL explicitly.', 'SEVERE')
                return False
        elif type(calibratorServiceURL)!=str or len(calibratorServiceURL)==0:
            casalog.post('Parameter calibratorServiceURL value must be None or a non-empty string, e.g. "https://almascience.org/sc/flux"', 'SEVERE')
            return False
            
        print("Will use calibratorServiceURL = '"+str(calibratorServiceURL)+"' in calls to aU.getALMAFluxForMS().")


    mytb = aU.createCasaTool(tbtool)
    mymsmd = msmdtool()


    ######################################
    # parse input parameters

    # useLocalAlmaHelper
    if (useLocalAlmaHelper):
        try:
            from almahelpers_localcopy import tsysspwmap2
            casalog.post('Using tsysspwmap2() from almahelpers_localcopy for generating Tsys SPW maps!', 'INFO')
            print('*** Using tsysspwmap2() from almahelpers_localcopy for generating Tsys SPW maps! ***')
        except:
            casalog.post('Module almahelpers_localcopy is not available. Please set useLocalAlmaHelper=False .', 'SEVERE')
            return False
    else:
        if mycasaversion < '5.9.9':
            from recipes.almahelpers import tsysspwmap
        else:
            from casarecipes.almahelpers import tsysspwmap

    # step
    availableSteps = ['calib', 
                      'fluxcal', 
                      'wvr', 
                      'calsurvey', 
                      'SDeff', 
                      'SDcalibLine', 
                      'SDcalibCont', 
                      'SDscience', 
                      'SDampcal']
    if step not in availableSteps:
        casalog.post("Step "+str(step)+" not valid.  Available values: "+str(availableSteps), 'SEVERE')
        return False

    if step in ['SDcalibLine', 'SDcalibCont', 'SDampcal', 'SDscience']:
        with_pointing_correction = True
    else:
        with_pointing_correction = False

    # refant
    if (type(refant) != str):
        casalog.post("refant must be a string", 'SEVERE')
        return False

    # splitMyScienceSpw, reindexMyScienceSpw
    if splitMyScienceSpw and (not reindexMyScienceSpw) and mycasaversion < '5.4':
        casalog.post("splitMyScienceSpw = True and reindexMyScienceSpw = False is not supported in CASA versions < 5.4", 'SEVERE')
        return False

    # state
    if state != '':
        casalog.post('The parameter "state" is no longer supported. Please set it to empty string.', 'SEVERE')
        return False

    # queue
    if queue != '':
        casalog.post('The parameter "queue" is no longer supported. Please set it to empty string.', 'SEVERE')
        return False

    # msNames
    if type(msNames) == str: 
        msNames = [msNames]

    # projectCode, schedblockName, schedblockUid, step
    if ((projectCode != '' and schedblockName != '') or schedblockUid != '') and step in ['calib', 'SDscience']:
        casalog.post('Automatic EB determination based on projectCode and schedblockName or schedblockUid is no longer supported.\n'\
                     +'Please set these parameters to empty string.', 'SEVERE')
        return False

    if type(allowHybrid) != bool:
        casalog.post('Parameter allowHybrid must be True or False.', 'SEVERE')
        return False

    if type(combineB2BLFspws) != bool:
        casalog.post('Parameter combineB2BLFspws must be True or False.', 'SEVERE')
        return False

    if type(combineB2BLFHFspws) != bool:
        casalog.post('Parameter combineB2BLFHFspws must be True or False.', 'SEVERE')
        return False

    if type(combineB2BDGCspws) != bool:
        casalog.post('Parameter combineB2BDGCspws must be True or False.', 'SEVERE')
        return False

    if type(includeRenorm) != bool:
        casalog.post('Parameter includeRenorm must be True or False.', 'SEVERE')
        return False

    if type(legacyfixes) != bool:
        casalog.post('Parameter legacyfixes must be True or False.', 'SEVERE')
        return False

    phaseDiff=False # may be overridden later in the case of BWSW or B2B

    # Logic for triggering various sections for B2B CASEs
    # Case A - Vanilla no combine
    # Case B - LF combine only in gaincal - combineB2BLFspws
    # Case C - LF and HF combine, in bandpass, and other gaincal for phase solns (_not_ B2B offset)
    #         - activate also the other mode as it is an 'add-on' 
    if combineB2BLFHFspws:
        combineB2BLFspws = True

    # Case D - LF and HF combine PLUS also combine the Diffgain cal solutions in the B2B offset
    #      - activate all other cases as this is an 'add-on'
    if combineB2BDGCspws:
        combineB2BLFHFspws = True
        combineB2BLFspws = True

    # With this CASE logic sequence, only the minimal 'trigger' needs to be coded in, i.e. combineB2BLF only where required
    # and will trigger in all other cases.
    
    # end initial parameter parsing
    ##################################

    currDir = os.getcwd()

    asis1 = 'Antenna Station Receiver Source CalAtmosphere CalWVR CorrelatorMode SBSummary'

    asdmNames = []
    for msName in msNames:
        myname = msName.rstrip(os.path.sep)
        if re.search('^uid\:\/\/[0-9a-z]+\/[0-9a-z]+\/[0-9a-z]+$', myname, re.IGNORECASE) is not None: # this is an ASDM name with slashes
            myname = re.sub(':|/', '_', myname) # convert to underscore notation

        if myname[-3:] == '.ms': # this has an ms extension    
            myname = myname[0:-3] # remove it
        
        if re.search('^uid___[0-9a-z]+_[0-9a-z]+_[0-9a-z]+$', myname, re.IGNORECASE) is not None: # this is now an ASDM name with underscores
            if not os.path.exists(myname) and not os.path.exists(myname+'.ms') and not os.path.exists('../'+myname):
                casalog.post('ERROR: '+myname+': Neither the asdm nor the ms exists. Please retrieve the ASDM before starting generateReducScript.', 'SEVERE')
                return False
            else:
                asdmNames.append(myname)
        else:
            casalog.post('ERROR: '+msName+' does not follow the ASDM naming convention.', 'SEVERE')
            return False

    print('Identified ASDMs '+str(asdmNames))

    msNames=[]
    for asdmName in asdmNames: # create the MS for each ASDM if it doesn't exist, yet

        msName = asdmName+'.ms'
        msNames.append(msName)

        if not os.path.exists(msName): # create the MS

            if not os.path.exists(asdmName):
                asdmName='../'+asdmName
                if not os.path.exists(asdmName):
                    casalog.post('ERROR: '+msName+': Neither the asdm nor the ms exists. Please retrieve the ASDM before starting generateReducScript.', 'SEVERE')
                    return False

            print('NOTE: The asdm '+asdmName+' exists, but the ms does not exist, running importasdm.')

            importasdm(asdmName, vis=msName, asis=asis1, bdfflags=bdfflags, lazy=lazy, process_caldevice=False, with_pointing_correction=with_pointing_correction)

        print('Determining the date of the observation')
        mymsmd.open(msName)
        thestartmjd = mymsmd.timerangeforobs(0)['begin']['m0']['value']
        mymsmd.close()
        if thestartmjd < 57296 and not legacyfixes:
            legacyfixes = True
            print('** Data was taken on MJD '+str(round(thestartmjd,1))+', i.e. before 1 Oct 2015. Will generate code for legacyfixes.')
            
        
    msNames = sorted(msNames)
    valueMaps = {}

    for msName in msNames:

        if legacyfixes:
            if step not in ['fluxcal']: sfsdr.fixForCSV2555(msName)

        if os.getlogin == 'aod':
            sfsdr.listOfIntentsWithSources(msName)
            sfsdr.listobs3(msName, figfile=msName+'.listobs3.png')
            sfsdr.plotAntennas(msName)

        spwInfo = sfsdr.getSpwInfo(msName, caching=True)
        spwIds = sorted(spwInfo.keys())
        print('Value mapping for MS '+msName+' ...')
        vm = aU.ValueMapping(msName)
        valueMaps[msName] = vm
        spwScans = vm.getScansForSpw(spwIds[0]).tolist()
        for j in spwIds:
            if vm.getScansForSpw(j).tolist() != spwScans:
                print('WARNING: The scans are not the same for all science spws.')
                print('WARNING: The script generator is not compatible with this, it will very likely fail.')
                print('WARNING: If it does not fail, do not expect the reduction script to be good. Please check it carefully.')


        ### determine if B2B or BWSW dataset or none of the above

        isB2B, isBWSW = isB2BorBWSW(msName, valueMaps)


    ##### step 'wvr' #####
    if step == 'wvr':

        for msName in msNames:

            print("\n*** Working on "+msName+" **********************************************")

            f1 = open(msName+'.scriptForWVRCalibration.py', 'w')
            print("import re\n", file=f1)
            print("es = aU.stuffForScienceDataReduction() \n\n", file=f1)
            if mycasaversion < '5.1':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', casaVersion) == None:", file=f1)
            elif mycasaversion < '5.9':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in cu.version().tolist()[:-1]])) == None:", file=f1)
            else:
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in casalith.version()[:-1]])) == None:", file=f1)
            print(" sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: "+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"')\n\n", file=f1)

            myRefAnt = refant
            if myRefAnt == '': myRefAnt = sfsdr.getRefAntenna(msName)
            print("# Using reference antenna = "+myRefAnt+"\n", file=f1)

            print(doAprioriFlagging(msName, valueMaps=valueMaps), file=f1)
            print(doGenerateWVRCalTable(msName, valueMaps=valueMaps), file=f1)
            print("es.wvr_stat(ms1='"+msName+"', refAnt='"+myRefAnt+"', qa2_output_dir='./')\n", file=f1)

            f1.close()

    ##### step 'calib' #####
    if step in ['calib']:

        for msName in msNames:

            print("\n*** Working on "+msName+" **********************************************")

            mystepdict = {}
            mystepindent = "  "
            myrepgenpardict = {} # parameters for the final QA2 report generation are gathered here

            tsysmap = ''
            if re.search('^3.3', mycasaversion) == None and skipSyscalChecks == False: # do Syscall checks

                print("\n*** ANALYSIS OF TSYS TABLE ***")

                print("\n*** SEARCH FOR NEGATIVE TSYS ***")
                aU.detectNegativeTsys(vis = msName, edge = 8, showfield = True)

                print("\n*** SEARCH FOR NEGATIVE TREC ***")
                aU.detectNegativeTrx(vis = msName, edge = 8, showfield = True)

                os.system('rm -Rf '+msName+'.tsys.temp')
                gencal(vis = msName, caltable = msName+'.tsys.temp', caltype = 'tsys')

                print("\n*** SEARCH FOR MISSING SCANS IN SYSCAL TABLE ***")
                mymsmd.open(msName)

                scans1 = sorted(mymsmd.scansforintent('CALIBRATE_ATMOSPHERE#*').tolist()) # T. Hunter 2014-08-14
                # There are cases, e.g. uid___A002_X7ea111_Xc03.ms where the OFF_SOURCE intent is present
                # in the ms but the ON_SOURCE is not, but the corresponding Tsys value is present.
                # As of 4.2.0, the forintent methods of msmd accept the wildcard character.

                mymsmd.close()

                mytb.open(msName+'.tsys.temp')                    
                scans2 = sorted(np.unique(mytb.getcol('SCAN_NUMBER')).tolist())
                mytb.close()

                if len(scans1) == 0 or len(scans2) == 0 or scans1 != scans2:
                    print("len(scans1)=%d, len(scans2)=%d" % (len(scans1), len(scans2)))
                    casalog.post('ERROR: THE SYSCAL TABLE IS MISSING ONE (OR MORE) SCAN(S). IT MAY BE NECESSARY TO RE-GENERATE IT.', 'SEVERE')
                    return False
                else:
                    print("-> OK")

                if useLocalAlmaHelper == True:
                    tsysmap = tsysspwmap2(vis = msName, tsystable = msName+'.tsys.temp', tsysChanTol=tsysChanTol)

                os.system('rm -Rf '+msName+'.tsys.temp')

            ##############################
            # Start the calibration script
            f1 = open(msName+'.scriptForCalibration.py', 'w')

            print("import re\n", file=f1)
            print("import os\n", file=f1)
            if (mycasaversion < '5.9'):
                print("import casadef\n", file=f1)
            else:
                print("import casalith\n", file=f1)

            addRenormStep=includeRenorm

            if addRenormStep:
                if (mycasaversion < '6.2.1'):
                    casalog.post('You are running a CASA version earlier than 6.2.1. Renorm cannot be applied!', 'WARN') 
                    sys.exit('Use CASA >= 6.2.1 or set includeRenorm=False !')
                print('\n*** includeRenorm == True: preparing Renormalization step by importing ACreNorm module from pipeline ...')
                try:
                    from pipeline.extern.almarenorm import ACreNorm
                except:
                    casalog.post('Could not import ACreNorm from pipeline.extern.almarenorm !', 'WARN')
                    casalog.post('For ALMA QA2, you need to run with a CASA version containing the pipeline!', 'WARN')
                    casalog.post('For other purposes, you may run with parameter includeRenorm=False .', 'WARN')
                    sys.exit('Could not import ACreNorm from pipeline.extern.almarenorm ! Use CASA >= 6.2.1 with pipeline or set includeRenorm=False !')


            print("if applyonly != True: es = aU.stuffForScienceDataReduction() \n\n", file=f1)
            if mycasaversion < '5.1':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', casaVersion) is None:", file=f1)
            elif mycasaversion < '5.9':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in cu.version().tolist()[:-1]])) == None:", file=f1)
            else:
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in casalith.version()[:-1]])) == None:", file=f1)
            print(" sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: "+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"')\n\n", file=f1)

            # add the list of intents as comments
            print('Generating list of intents and fields ...')
            print(listOfIntentsWithFields(msName), file=f1)

            # determine reference antenna
            myRefAnt = refant
            if myRefAnt == '': myRefAnt = sfsdr.getRefAntenna(msName)
            print("\n# Using reference antenna = "+myRefAnt, file=f1)

            isFullP = isFullPol(msName, valueMaps)

            if isFullP:
                print("NOTE: This is a full polarisation dataset.")
                print("#     (will use it strictly since this is a full polarisation dataset)\n", file=f1)
            else:
                print("", file=f1)
                

            ### add importasdm step
            stext = "if os.path.exists('"+msName+"') == False:\n"

            stext += "  importasdm('"+re.findall('.+(?=\.ms)', msName, re.IGNORECASE)[0]+"', asis='"+asis1+"', bdfflags="+str(bdfflags)+", lazy="+str(lazy)+", process_caldevice=False)"

            stext += "\n  if not os.path.exists('"+msName+".flagversions'):\n    print('ERROR in importasdm. Output MS is probably not useful. Will stop here.')\n    thesteps = []"

            if legacyfixes:
                stext += "\nif applyonly != True: es.fixForCSV2555('"+msName+"')"
                
            sfsdr.addReducScriptStep(f1, mystepdict, "Import of the ASDM", stext, mystepindent, applyonly=True)

            if legacyfixes:
                ### add fixsyscaltimes step: see CAS-4981, CSV-2841, ICT-3642; sometimes needed for cycle 0-2.
                if mycasaversion < '5.9.9':
                    stext = "from recipes.almahelpers import fixsyscaltimes\nfixsyscaltimes(vis = '"+msName+"')"
                else:
                    stext = "from casarecipes.almahelpers import fixsyscaltimes\nfixsyscaltimes(vis = '"+msName+"')"
                sfsdr.addReducScriptStep(f1, mystepdict, "Fix of SYSCAL table times", stext, mystepindent, applyonly=True)

            print('print("# A priori calibration")\n', file=f1)

            if legacyfixes:
                ### add fixplanets step if needed 
                stext = doRunFixPlanets(msName)
                if stext is not None: 
                    sfsdr.addReducScriptStep(f1, mystepdict, "Running fixplanets on fields with 0,0 coordinates", stext, mystepindent, applyonly=True)
                else:
                    print('No (0,0) coordinates found. fixplanets step not needed.')

            ### add listobs step
            stext = "os.system('rm -rf %s.listobs')\n" %(msName) # Added by CLB
            stext += "listobs(vis = '"+msName+"',\n  listfile = '"+msName+".listobs')\n\n" # Modified by CLB
            sfsdr.addReducScriptStep(f1, mystepdict, "listobs", stext, mystepindent)

            ### add apriori flagging step
            stext = doAprioriFlagging(msName, valueMaps=valueMaps)
            sfsdr.addReducScriptStep(f1, mystepdict, "A priori flagging", stext, mystepindent)

            ### add wvrgcal step
            wvrCalTableName = []
            stext = doGenerateWVRCalTable(msName, wvrCalTableName, refant=myRefAnt, remcloud=remcloud, valueMaps=valueMaps, isB2B=isB2B, isBWSW=isBWSW)
            sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the WVR cal table", stext, mystepindent)

            ### add Tsys step
            tsysCalTableName = []
            stext = doGenerateTsysCalTable(msName, tsysCalTableName, isB2B=isB2B)
            sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the Tsys cal table", stext, mystepindent)

            ### add antpos (if needed) and first applycal step
            if corrAntPos == True:
                stext = sfsdr.correctMyAntennaPositions(msName)  ## use default setting for maxSearchDays, see also SCIREQ-2684
                if stext is not None:
                    sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the antenna position cal table", stext, mystepindent)
                    stext = doApplyAprioriCalTables(msName, tsys=tsysCalTableName[0], wvr=wvrCalTableName[0], antpos=msName+'.antpos', tsysmap=tsysmap, tsysChanTol=tsysChanTol, tsysPerField=tsysPerField, valueMaps=valueMaps)
                    sfsdr.addReducScriptStep(f1, mystepdict, "Application of the WVR, Tsys and antpos cal tables", stext, mystepindent, applyonly=True)
                else:
                    stext = doApplyAprioriCalTables(msName, tsys=tsysCalTableName[0], wvr=wvrCalTableName[0], tsysmap=tsysmap, tsysChanTol=tsysChanTol, tsysPerField=tsysPerField, valueMaps=valueMaps)
                    sfsdr.addReducScriptStep(f1, mystepdict, "Application of the WVR and Tsys cal tables", stext, mystepindent, applyonly=True)
            else:
                stext = doApplyAprioriCalTables(msName, tsys=tsysCalTableName[0], wvr=wvrCalTableName[0], tsysmap=tsysmap, tsysChanTol=tsysChanTol, tsysPerField=tsysPerField, valueMaps=valueMaps)
                sfsdr.addReducScriptStep(f1, mystepdict, "Application of the WVR and Tsys cal tables", stext, mystepindent, applyonly=True)

            ### add splitout step
            stext = doSplitOut(msName, splitMyScienceSpw=splitMyScienceSpw, timebin=timeBinForFinalData, reindexMyScienceSpw=reindexMyScienceSpw)
            mysteptitleaddon = ""
            if timeBinForFinalData > 0.:
                mysteptitleaddon = " and time average"
            sfsdr.addReducScriptStep(f1, mystepdict, "Split out science SPWs"+mysteptitleaddon, stext, mystepindent, applyonly=True)


            #### end of apriori calibration steps generation ####

            print('print("# Calibration")\n', file=f1)

            ### add listobs step
            stext = "os.system('rm -rf %s.split.listobs')\n" % (msName) # Added by CLB
            stext += "listobs(vis = '"+msName+".split',\n  listfile = '"+msName+".split.listobs')\n\n"
            stext += doSaveFlags(msName+'.split', name='Original')
            sfsdr.addReducScriptStep(f1, mystepdict, "Listobs, and save original flags", stext, mystepindent)

            ### add initial flagging step
            stext = doInitialFlagging(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw)
            sfsdr.addReducScriptStep(f1, mystepdict, "Initial flagging", stext, mystepindent)

            ### add setjy step
            stext, theFluxCalNames = doRunSetjy(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, 
                                                useCalibratorService=useCalibratorService, calibratorServiceURL=calibratorServiceURL,  
                                                isB2B=isB2B, valueMaps=valueMaps)
            sfsdr.addReducScriptStep(f1, mystepdict, "Putting a model for the flux calibrator(s)", stext, mystepindent)

            ### add bandpass step
            thebpassCalTableName = [bpassCalTableName]            
            if bpassCalTableName == '':
                stext = doSaveFlags(msName+'.split', name='BeforeBandpassCalibration')
                sfsdr.addReducScriptStep(f1, mystepdict, "Save flags before bandpass cal", stext, mystepindent, applyonly=True)
                thebpassCalTableName = []
                stext = doBandpassCalibration(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, 
                                              refant=myRefAnt, calTableName=thebpassCalTableName, lowSNR=lowSNR, lbc=lbc, 
                                              phaseDiff=phaseDiff, isB2B=isB2B, isBWSW=isBWSW, combineB2BLFHFspws=combineB2BLFHFspws, isFullP=isFullP, 
                                              bpassCalId=bpassCalId, theFluxCalNames=theFluxCalNames, vetoBPBootstrap=vetoBPBootstrap,
                                              valueMaps=valueMaps)
                #### NOTE: thebpassCalTableName[0] is set by doBandpassCalibration!
                sfsdr.addReducScriptStep(f1, mystepdict, "Bandpass calibration", stext, mystepindent)

            ### add saveflags step
            stext = doSaveFlags(msName+'.split', name='BeforeGainCalibration')
            sfsdr.addReducScriptStep(f1, mystepdict, "Save flags before gain cal", stext, mystepindent, applyonly=True)

            ### add gain calibration step
            phaseDiffCalTableName = []
            ampForSci = []
            if isB2B:
                stext = doB2BGainCalibrationPartI(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, refant=myRefAnt, 
                                                  bandpass=thebpassCalTableName[0], valueMaps=valueMaps, combineB2BLFspws=combineB2BLFspws,
                                                  combineB2BLFHFspws=combineB2BLFHFspws, combineB2BDGCspws=combineB2BDGCspws)
                sfsdr.addReducScriptStep(f1, mystepdict, "B2B Gain calibration Part I", stext, mystepindent)

                stext = doB2BGainCalibrationPartII(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, refant=myRefAnt, 
                                                   bandpass=thebpassCalTableName[0], ampForSci=ampForSci, 
                                                   valueMaps=valueMaps, combineB2BLFspws=combineB2BLFspws,
                                                   combineB2BLFHFspws=combineB2BLFHFspws, combineB2BDGCspws=combineB2BDGCspws)
                sfsdr.addReducScriptStep(f1, mystepdict, "B2B Gain calibration Part II", stext, mystepindent)

            else:
                stext = doGainCalibration(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, refant=myRefAnt, 
                                          bandpass=thebpassCalTableName[0], phaseDiffCalTableName=phaseDiffCalTableName, ampForSci=ampForSci, 
                                          phaseDiff=phaseDiff, 
                                          isBWSW=isBWSW,
                                          valueMaps=valueMaps)
                sfsdr.addReducScriptStep(f1, mystepdict, "Gain calibration", stext, mystepindent)

            ## for CASA 6.6.1 and later: add renorm step before applycal
            renormTableName=''
            renorm_message=''
            if addRenormStep and (mycasaversion > '6.6.0'):
                # we have already checked above that the present CASA version actually has the renorm recipes available
                try:
                    rn = ACreNorm(msName)
                except:
                    print(sys.exc_info())
                    casalog.post('Error in ACreNorm from pipeline.extern.almarenorm ! You may want to run again with includeRenorm=False.', 'WARN')
                    sys.exit('Error in ACreNorm from pipeline.extern.almarenorm ! You may want to run again with includeRenorm=False.')
                if rn.tdm_only:
                    casalog.post('This dataset only has TDM SPWs. No renorm investigation needed.', 'INFO')
                else:
                    stext, renormTableName = doRenormTable(msName, msName1=msName+'.split', isB2B=isB2B, isBWSW=isBWSW,
                                                           iHaveSplitMyScienceSpw=reindexMyScienceSpw, valueMaps=valueMaps)
                    renorm_message = ' and the renorm table'
                    sfsdr.addReducScriptStep(f1, mystepdict, "Renormalization", stext, mystepindent)


            ### add safeflags step
            stext = doSaveFlags(msName+'.split', name='BeforeApplycal')
            sfsdr.addReducScriptStep(f1, mystepdict, "Save flags before applycal", stext, mystepindent, applyonly=True)

              
            ### add final applycal step
            if isB2B:
                stext = doApplyB2BBandpassAndGainCalTables(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, 
                                                           bandpass=thebpassCalTableName[0], valueMaps=valueMaps, combineB2BLFspws=combineB2BLFspws,
                                                           combineB2BLFHFspws=combineB2BLFHFspws, combineB2BDGCspws=combineB2BDGCspws,
                                                           renorm=renormTableName)
            else:
                stext = doApplyBandpassAndGainCalTables(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, 
                                                        bandpass=thebpassCalTableName[0], phaseForCal=msName+'.split.phase_int', 
                                                        phaseForSci=msName+'.split.phase_inf', flux=msName+'.split.flux_inf', 
                                                        phaseDiffCalTableName=phaseDiffCalTableName, ampForSci=ampForSci, 
                                                        valueMaps=valueMaps, renorm=renormTableName)

            sfsdr.addReducScriptStep(f1, mystepdict, "Application of the bandpass and gain cal tables"+renorm_message, stext, mystepindent, applyonly=True)


            if addRenormStep and (mycasaversion < '6.6.1'): ### add a renormalisation step to investigate and fix problems with ATM lines, if there are FDM SPWs
                # we have already checked above that the present CASA version actually has the renorm recipes available
                try:
                    rn = ACreNorm(msName)
                except:
                    print(sys.exc_info())
                    casalog.post('Error in ACreNorm from pipeline.extern.almarenorm ! You may want to run again with includeRenorm=False.', 'WARN')
                    sys.exit('Error in ACreNorm from pipeline.extern.almarenorm ! You may want to run again with includeRenorm=False.')
                if rn.tdm_only:
                    casalog.post('This dataset only has TDM SPWs. No renorm investigation needed.', 'INFO')
                else:
                    stext = doRenorm(msName, msName1=msName+'.split', isB2B=isB2B, isBWSW=isBWSW,
                                     iHaveSplitMyScienceSpw=reindexMyScienceSpw, valueMaps=valueMaps)
                    sfsdr.addReducScriptStep(f1, mystepdict, "Run renormalization", stext, mystepindent, applyonly=True)

            ### add final splitout step 
            stext = doSplitOut(msName, msName1=msName+'.split', outMsName=msName+'.split.cal', allowHybrid=allowHybrid, intentsToDiscard='ATMOSPHERE|POINTING', iHaveSplitMyScienceSpw=reindexMyScienceSpw)
            sfsdr.addReducScriptStep(f1, mystepdict, "Split out corrected column", stext, mystepindent, applyonly=True)

            ### add final flagsave step
            stext = doSaveFlags(msName+'.split.cal', name='AfterApplycal')
            sfsdr.addReducScriptStep(f1, mystepdict, "Save flags after applycal", stext, mystepindent, applyonly=True)

            ### add QA2 report generation step

            stext = doQa2ReportGeneration(msName, refant=myRefAnt, isB2B=isB2B, isBWSW=isBWSW, iHaveSplitMyScienceSpw=reindexMyScienceSpw)
            sfsdr.addReducScriptStep(f1, mystepdict, "Generate QA2 Report", stext, mystepindent, applyonly=True)


            ### finish script by adding header
            sfsdr.prependReducScriptHeader(f1, mystepdict, "Created using "+version(), mystepindent)

            f1.close()


    ########################################

    if step == 'fluxcal':

        if os.path.exists('allFluxes.txt') == False:
            sfsdr.generateFluxFile(msNames)
        else:
            print('File allFluxes.txt already exists, it will be loaded.')

        myRefAnt = refant # note: no need to run getRefAntenna because it will be called in doFluxCalibration if refant=''
        f1 = open('scriptForFluxCalibration.py', 'w')
        print("import re\n", file=f1)
        if mycasaversion < '5.1':
            print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', casaVersion) == None:", file=f1)
        elif mycasaversion < '5.9':
            print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in cu.version().tolist()[:-1]])) == None:", file=f1)
        else:
            print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in casalith.version()[:-1]])) == None:", file=f1)
        print(" sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: "+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"')\n\n", file=f1)
        print(doFluxCalibration(msNames, refant=myRefAnt), file=f1)
        f1.close()

    ######################################


    if step == 'calsurvey':
        mystepdict = {}
        mystepindent = "  "

        for msName in msNames:

            print("\n*** Working on "+msName+" **********************************************")

            tsysmap = ''
            if re.search('^3.3', mycasaversion) == None and skipSyscalChecks == False:

                print("\n*** ANALYSIS OF TSYS TABLE ***")

                print("\n*** SEARCH FOR NEGATIVE TSYS ***")
                aU.detectNegativeTsys(vis = msName, edge = 8, showfield = True)

                print("\n*** SEARCH FOR NEGATIVE TREC ***")
                aU.detectNegativeTrx(vis = msName, edge = 8, showfield = True)

                os.system('rm -Rf '+msName+'.tsys.temp')
                gencal(vis = msName, caltable = msName+'.tsys.temp', caltype = 'tsys')

                print("\n*** SEARCH FOR MISSING SCANS IN SYSCAL TABLE ***")

                mymsmd = msmdtool()
                mymsmd.open(msName)

                scans1 = mymsmd.scansforintent('CALIBRATE_ATMOSPHERE*')

                mymsmd.close()

                mytb.open(msName+'.tsys.temp')                    
                scans2 = np.unique(mytb.getcol('SCAN_NUMBER'))
                mytb.close()

                if (scans1 == scans2).all():
                    print("-> OK")
                else:
                    casalog.post('ERROR: THE SYSCAL TABLE IS MISSING ONE (OR MORE) SCAN(S). IT MAY BE NECESSARY TO RE-GENERATE IT.', 'SEVERE')
                    return False

                if useLocalAlmaHelper == True:
                    tsysmap = tsysspwmap2(vis = msName, tsystable = msName+'.tsys.temp', tsysChanTol=tsysChanTol)

                os.system('rm -Rf '+msName+'.tsys.temp')

            f1 = open(msName+'.scriptForCalibration.py', 'w')
            print("import re\n", file=f1)
            print("es = aU.stuffForScienceDataReduction() \n\n", file=f1)
            if mycasaversion < '5.1':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', casaVersion) == None:", file=f1)
            elif mycasaversion < '5.9':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in cu.version().tolist()[:-1]])) == None:", file=f1)
            else:
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in casalith.version()[:-1]])) == None:", file=f1)
            print(" sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: "+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"')\n\n", file=f1)

            myRefAnt = refant
            if myRefAnt == '': myRefAnt = sfsdr.getRefAntenna(msName)
            print("# Using reference antenna = "+myRefAnt+"\n", file=f1)

            print('print("# A priori calibration")\n', file=f1)
            if legacyfixes:
                stext = doRunFixPlanets(msName)
                if stext is not None: sfsdr.addReducScriptStep(f1, mystepdict, "Running fixplanets on fields with 0,0 coordinates", stext, mystepindent)

            stext = "os.system('rm -rf %s.listobs')\n" %(msName) # Added by CLB
            stext += "listobs(vis = '"+msName+"',\n  listfile = '"+msName+".listobs')\n\n" # Modified by CLB
            sfsdr.addReducScriptStep(f1, mystepdict, "listobs", stext, mystepindent)
            stext = doAprioriFlagging(msName, valueMaps=valueMaps)
            sfsdr.addReducScriptStep(f1, mystepdict, "A priori flagging", stext, mystepindent)
            wvrCalTableName = []
            stext = doGenerateWVRCalTable(msName, wvrCalTableName, valueMaps=valueMaps, isB2B=isB2B)
            sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the WVR cal table", stext, mystepindent)
            tsysCalTableName = []
            stext = doGenerateTsysCalTable(msName, tsysCalTableName, isB2B=isB2B)
            sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the Tsys cal table", stext, mystepindent)

            if corrAntPos == True:
                stext = sfsdr.correctMyAntennaPositions(msName) # use default setting for maxSearchDays, see also SCIREQ-2684
            if corrAntPos == True and stext is not None:
                sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the antenna position cal table", stext, mystepindent)
                stext = doApplyAprioriCalTables(msName, tsys=tsysCalTableName[0], wvr=wvrCalTableName[0], antpos=msName+'.antpos', tsysmap=tsysmap, valueMaps=valueMaps)
                sfsdr.addReducScriptStep(f1, mystepdict, "Application of the WVR, Tsys and antpos cal tables", stext, mystepindent)
            else:
                stext = doApplyAprioriCalTables(msName, tsys=tsysCalTableName[0], wvr=wvrCalTableName[0], valueMaps=valueMaps)
                sfsdr.addReducScriptStep(f1, mystepdict, "Application of the WVR and Tsys cal tables", stext, mystepindent)

            stext = doSplitOut(msName, splitMyScienceSpw=splitMyScienceSpw, timebin=timeBinForFinalData, reindexMyScienceSpw=reindexMyScienceSpw) 
            mysteptitleaddon = ""
            if timeBinForFinalData > 0.:
                mysteptitleaddon = " and time average"
            sfsdr.addReducScriptStep(f1, mystepdict, "Split out science SPWs"+mysteptitleaddon, stext, mystepindent)

            print('print("# Calibration")\n', file=f1)
            stext = "os.system('rm -rf %s.split.listobs')\n" % (msName)
            stext += "listobs(vis = '"+msName+".split',\n  listfile = '"+msName+".split.listobs')\n\n" \
                + doClearPointingTable(msName+'.split') \
                + doSaveFlags(msName+'.split', name='Original')
            sfsdr.addReducScriptStep(f1, mystepdict, "Listobs, clear pointing table, and save original flags", stext, mystepindent)
            stext = doInitialFlagging(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw)
            sfsdr.addReducScriptStep(f1, mystepdict, "Initial flagging", stext, mystepindent)
            stext, theFluxCalNames = doRunSetjy(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, 
                                                useCalibratorService=useCalibratorService, calibratorServiceURL=calibratorServiceURL, valueMaps=valueMaps)
            sfsdr.addReducScriptStep(f1, mystepdict, "Putting a model for the flux calibrator(s)", stext, mystepindent)
            stext = doSaveFlags(msName+'.split', name='BeforeBandpassCalibration')
            sfsdr.addReducScriptStep(f1, mystepdict, "Save flags before bandpass cal", stext, mystepindent)
            thebpassCalTableName = []
            stext = doBandpassCalibration(msName, msName1=msName+'.split', bpassCalId=bpassCalId, iHaveSplitMyScienceSpw=reindexMyScienceSpw, refant=myRefAnt, 
                                          calTableName=thebpassCalTableName, theFluxCalNames=theFluxCalNames, vetoBPBootstrap=vetoBPBootstrap,
                                          valueMaps=valueMaps)
            sfsdr.addReducScriptStep(f1, mystepdict, "Bandpass calibration", stext, mystepindent)
            stext = doSaveFlags(msName+'.split', name='BeforeGainCalibration')
            sfsdr.addReducScriptStep(f1, mystepdict, "Save flags before gain cal", stext, mystepindent)
            stext = doGainCalibration(msName, msName1=msName+'.split', iHaveSplitMyScienceSpw=reindexMyScienceSpw, refant=myRefAnt,
                                      bandpass=thebpassCalTableName[0], gaintypeForAmp='T', valueMaps=valueMaps)
            sfsdr.addReducScriptStep(f1, mystepdict, "Gain calibration", stext, mystepindent)

            sfsdr.prependReducScriptHeader(f1, mystepdict, "Created using "+version(), mystepindent)

            f1.close()


    if step == 'SDeff':

        for msName in msNames:

            print("\n*** Working on "+msName+" **********************************************")

            fieldNames = sfsdr.getIntentsAndSourceNames(msName)['OBSERVE_TARGET']['name']
            fieldNames = sorted(dict.fromkeys(fieldNames).keys())
            if len(fieldNames) != 1: 
                casalog.post('ERROR: Unexpected number of fields.', 'SEVERE')
                return False
            fieldNames = fieldNames[0]
            if fieldNames.upper() in ['VENUS', 'MARS', 'JUPITER', 'URANUS', 'NEPTUNE', 'IO', 'EUROPA', 'GANYMEDE', 'CALLISTO', 'TITAN', 'CERES', 'JUNO', 'PALLAS', 'VESTA', 'HYGEIA']:
                print('Observation type = SSO')
                sdEffType = 'SSO'
            else:
                print('Observation type = QSO')
                sdEffType = 'QSO'

            mytb.open(msName+'/OBSERVATION')
            obsTimeRange = mytb.getcol('TIME_RANGE')
            obsTime = (obsTimeRange[0]+obsTimeRange[1])/2.0
            obsTime = ((obsTime/86400.0)+2400000.5-2440587.5)*86400.0
            obsTime = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTime[0])))
            mytb.close()

            mystepdict = {}
            mystepindent = "  "

            f1 = open(msName+'.scriptForSDefficiencies.py', 'w')
            print("import re\n", file=f1)
            print("es = aU.stuffForScienceDataReduction()\n", file=f1)
            print("import analysisUtilsForSD as aUsd\n\n", file=f1)
            if mycasaversion < '5.1':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', casaVersion) == None:", file=f1)
            elif mycasaversion < '5.9':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in cu.version().tolist()[:-1]])) == None:", file=f1)
            else:
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in casalith.version()[:-1]])) == None:", file=f1)
            print(" sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: "+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"')\n\n", file=f1)

            asdmName1 = re.findall('^uid___[0-9a-z]+_[0-9a-z]+_[0-9a-z]+', msName, re.IGNORECASE)[0]

            stext = "aUsd.continuumReducer2('"+asdmName1+"')"
            sfsdr.addReducScriptStep(f1, mystepdict, "Make continuum images", stext, mystepindent)

            mytb.open(msName+'/ANTENNA')
            antNames = mytb.getcol('NAME')
            mytb.close()

            antNames1 = []
            for i in antNames:
                if re.search('^CM[0-9]+', i) == None: antNames1.append(i)
            if len(antNames1) == 0: 
                casalog.post('ERROR: No antenna to process.', 'SEVERE')
                return False
                
            if msName in valueMaps.keys():
                vm = valueMaps[msName]
                print('Using canned ValueMap.')
            else:
                vm = aU.ValueMapping(msName)
                valueMaps[msName] = vm

            sciScans = vm.getScansForIntent('OBSERVE_TARGET#ON_SOURCE').tolist()
            sciScans = [str(i) for i in sciScans]
            sciScans = ','.join(sciScans)

            weatherInfo = {}
            for i in antNames1:
                weatherInfo[i] = aU.getWeather(msName, antenna=i, scan=sciScans, getSolarDirection=False)[0]

            spwInfo = sfsdr.getSpwInfo(msName, caching=True)
            spwIds = sorted(spwInfo.keys())

            stext = ''

            for i in antNames1:
                for j in spwIds:

                    msName1 = asdmName1 + '.' + i + '.cal.ms'
                    imgName1 = asdmName1 + '.' + i + '.SPW' + str(j) + '.SF.im'

                    stext += "aUsd.SdImDWtApp(ms_SD = '"+msName1+"',\n  im_SD = '"+imgName1+"',\n  outimage = '"+imgName1+".wt',\n  spwid = "+str(j)+",\n  gridfunc = 'SF')\n"

                stext += '\n'

            sfsdr.addReducScriptStep(f1, mystepdict, "Mitigating noisy effect at edge of SD image", stext, mystepindent)

            if sdEffType == 'SSO':

                stext = 'ssoParams = {}\n\n'

                for i in antNames1:

                    msName1 = asdmName1 + '.' + i + '.cal.ms'

                    stext += "ssoParams['"+i+"'] = {}\n"
                    stext += "s = aUsd.sso_params('"+msName1+"')\n\n"

                    stext += "for i in "+str(spwIds)+":\n\n"
                    stext += "  s.setSpwId(int(i))\n"
                    stext += "  s.doCalc()\n"
                    stext += "  ssoResults = s.getResults()\n\n"
                    stext += "  ssoParams['"+i+"'][i] = {}\n"
                    stext += "  ssoParams['"+i+"'][i]['eqsize'] = ( ssoResults['Apparent Size'][0][0] + ssoResults['Apparent Size'][1][0] ) / 2.\n"
                    stext += "  ssoParams['"+i+"'][i]['psize'] = ( ssoResults['Apparent Size'][0][1] + ssoResults['Apparent Size'][1][1] ) / 2.\n"
                    stext += "  ssoParams['"+i+"'][i]['pnang'] = ( ssoResults['Apparent Size'][0][2] + ssoResults['Apparent Size'][1][2] ) / 2.\n"
                    stext += "  ssoParams['"+i+"'][i]['btemp'] = ( ssoResults['Brightness Temperature'][0] + ssoResults['Brightness Temperature'][1] ) / 2.\n\n"

                stext += "f = open('"+msName+".ssoParams.txt', 'w')\n"
                stext += "f.write(ssoParams+'\\n')\n"
                stext += "f.close()\n\n"

                sfsdr.addReducScriptStep(f1, mystepdict, "Obtain apparent size and brightness temperature for Solar system object", stext, mystepindent)

            stext = 'sdEffs = {}\n\n'

            for i in antNames1:

                msName1 = asdmName1 + '.' + i + '.cal.ms'

                stext += "for i in "+str(spwIds)+":\n\n"
                stext += "  a = aUsd.analysis_sdim('"+asdmName1+"."+i+".SPW'+str(i)+'.SF.im.wt')\n"

                if sdEffType == 'SSO':
                    stext += "  a.setSSOParams(eqsize = ssoParams['"+i+"'][i]['eqsize'],\n    psize = ssoParams['"+i+"'][i]['psize'],\n    pnang = ssoParams['"+i+"'][i]['pnang'],\n    btemp = ssoParams['"+i+"'][i]['btemp'])\n"
                    stext += "  a.doSSOAnalysis(antbeam = True)\n"
                    stext += "  results = a.getSSOResults()\n"
                    stext += "  #a.showSSOModel()\n"
                    stext += "  #a.showSSOResults()\n\n"
                else:
                    stext += "  a.setQSOFlux("+str(sdQSOflux)+")\n"
                    stext += "  a.doQSOAnalysis()\n"
                    stext += "  results = a.getQSOResults()\n"
                    stext += "  #a.showSdBeam()\n\n"

                stext += "  ij = len(sdEffs)\n"
                stext += "  sdEffs[ij] = {}\n"
                stext += "  sdEffs[ij]['execBlockUid'] = '"+asdmName1+"'\n"
                stext += "  sdEffs[ij]['obsTime'] = '"+obsTime+"'\n"
                stext += "  sdEffs[ij]['antennaName'] = '"+i+"'\n"
                stext += "  sdEffs[ij]['spwId'] = i\n"
                stext += "  sdEffs[ij]['frequency'] = results['Frequency']/1e9\n"
                stext += "  sdEffs[ij]['meanElevation'] = "+str(weatherInfo[i]['elevation'])+"\n"
                stext += "  sdEffs[ij]['meanTemp'] = "+str(weatherInfo[i]['temperature'])+"\n"
                stext += "  sdEffs[ij]['meanWindSpeed'] = "+str(weatherInfo[i]['windspeed'])+"\n"
                stext += "  sdEffs[ij]['effectiveBeamSize'] = results['Effective Beam Size'].tolist()\n"
                stext += "  sdEffs[ij]['mainBeamEfficiency'] = results['Main Beam Efficiency']\n\n"

            stext += "f = open('"+msName+".sdEfficiencies.txt', 'w')\n"
            stext += "f.write(sdEffs+'\\n')\n"
            stext += "f.close()\n\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Obtain efficiencies and (effective) beam size", stext, mystepindent)

            sfsdr.prependReducScriptHeader(f1, mystepdict, "Calculation of SD efficiencies\n# Created using "+version(), mystepindent)

            f1.close()

    if step in ['SDcalibLine', 'SDcalibCont', 'SDampcal', 'SDscience']:

        if step == 'SDscience' and len(msNames) > 1:

            spwInfo = sfsdr.getSpwInfo(msNames[0], caching=True)
            spwIds = sorted(spwInfo.keys())

            for j in range(1, len(msNames)):

                spwInfo1 = sfsdr.getSpwInfo(msNames[j], caching=True)
                spwIds1 = sorted(spwInfo1.keys())
                if spwIds1 != spwIds: 
                    print('WARNING: THE SCIENCE SPWS ARE NOT THE SAME FOR ALL EXECUTIONS.')

        imagingParams = {}

        for msName in msNames:

            print("\n*** Working on "+msName+" **********************************************")

            mystepdict = {}
            mystepindent = "  "

            tsysmap = ''
            if re.search('^3.3', mycasaversion) is None:

                print("\n*** ANALYSIS OF TSYS TABLE ***")

                print("\n*** SEARCH FOR NEGATIVE TSYS ***")
                aU.detectNegativeTsys(vis = msName, edge = 8, showfield = True)

                print("\n*** SEARCH FOR NEGATIVE TREC ***")
                aU.detectNegativeTrx(vis = msName, edge = 8, showfield = True)

                os.system('rm -Rf '+msName+'.tsys.temp')
                gencal(vis = msName, caltable = msName+'.tsys.temp', caltype = 'tsys')

                print("\n*** SEARCH FOR MISSING SCANS IN SYSCAL TABLE ***")
                mymsmd = msmdtool()
                mymsmd.open(msName)

                scans1 = mymsmd.scansforintent('CALIBRATE_ATMOSPHERE*')

                if ('CALIBRATE_ATMOSPHERE#REFERENCE' in mymsmd.intents()):
                    # The presence of extra AtmCals with a REFERENCE intent will cause scans1 != scans2
                    scansWithZeroLevel = mymsmd.scansforintent('CALIBRATE_ATMOSPHERE#REFERENCE')
                    scans1 = np.array(sorted(list(set(scans1)-set(scansWithZeroLevel))))
                mymsmd.close()

                mytb.open(msName+'.tsys.temp')                    
                scans2 = np.unique(mytb.getcol('SCAN_NUMBER'))
                mytb.close()

                if (scans1 == scans2).all():
                    print("-> OK")
                else:
                    print("scans1 = ", scans1)
                    print("scans2 = ", scans2)
                    casalog.post('ERROR: THE SYSCAL TABLE IS MISSING ONE (OR MORE) SCAN(S). IT MAY BE NECESSARY TO RE-GENERATE IT.', 'SEVERE')
                    return False

                if useLocalAlmaHelper == True:
                    tsysmap = tsysspwmap2(vis = msName, tsystable = msName+'.tsys.temp', tsysChanTol=tsysChanTol)

                os.system('rm -Rf '+msName+'.tsys.temp')

            if step == 'SDampcal':
                f1name = msName+'.scriptForSDampcalReduction.py'
            else:
                f1name = msName+'.scriptForSDCalibration.py'
            f1 = open(f1name, 'w')

            if(step != 'SDscience' or  mycasaversion < '6.4.3'): # add some functions for scaleAutocorr
                print("import os", file=f1)
                print("import re\n", file=f1)

                sourcelines1 = ''.join(inspect.getsourcelines(aU.createCasaTool)[0])
                if re.search('""".*?"""', sourcelines1, re.DOTALL) is not None:
                    sourcelines2 = re.findall('""".*?"""', sourcelines1, re.DOTALL)[0]
                    sourcelines1 = sourcelines1.replace(sourcelines2, '')
                    print(sourcelines1, file=f1)
                    
                sourcelines1 = ''.join(inspect.getsourcelines(aU.getDataColumnName)[0])
                if re.search('""".*?"""', sourcelines1, re.DOTALL) is not None:
                    sourcelines2 = re.findall('""".*?"""', sourcelines1, re.DOTALL)[0]
                    sourcelines1 = sourcelines1.replace(sourcelines2, '')
                    print(sourcelines1, file=f1)

                sourcelines1 = ''.join(inspect.getsourcelines(aU.scaleAutocorr)[0])
                if re.search('""".*?"""', sourcelines1, re.DOTALL) is not None:
                    sourcelines2 = re.findall('""".*?"""', sourcelines1, re.DOTALL)[0]
                    sourcelines1 = sourcelines1.replace(sourcelines2, '')
                    print(sourcelines1, file=f1)

            print("if applyonly != True: es = aU.stuffForScienceDataReduction()\n", file=f1)
            if mycasaversion < '5.1':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', casaVersion) == None:", file=f1)
            elif mycasaversion < '5.9':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in cu.version().tolist()[:-1]])) == None:", file=f1)
            else:
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in casalith.version()[:-1]])) == None:", file=f1)
            print("  sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: "+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"')\n\n", file=f1)

            print('Generating list of intents and fields ...')
            print(listOfIntentsWithFields(msName), file=f1)

            print("\n", file=f1)

            stext = "if os.path.exists('"+msName+"') == False:\n"

            stext += "  importasdm('"+re.findall('.+(?=\.ms)', msName, re.IGNORECASE)[0]+"', asis='"+asis1+"', bdfflags="+str(bdfflags)+", lazy="+str(lazy)+", process_caldevice=False, with_pointing_correction="+str(with_pointing_correction)+")"

            if legacyfixes:
                stext += "\nif applyonly != True: es.fixForCSV2555('"+msName+"')"
                
            sfsdr.addReducScriptStep(f1, mystepdict, "Import of the ASDM", stext, mystepindent, applyonly=True)

            if legacyfixes:
                stext = doRunFixPlanets(msName)
                if stext is not None: sfsdr.addReducScriptStep(f1, mystepdict, "Running fixplanets on fields with 0,0 coordinates", stext, mystepindent)

            mytb.open(msName+'/ANTENNA')
            antNames = mytb.getcol('NAME').tolist()
            mytb.close()

            stext = "os.system('rm -rf %s.listobs')\n" %(msName) # Added by CLB
            stext += "listobs(vis = '"+msName+"',\n  listfile = '"+msName+".listobs')\n\n" # Modified by CLB

            stext += "if applyonly != True:\n"
            stext += "  aU.getTPSampling(vis = '"+msName+"', showplot = True, plotfile = '"+msName+".sampling.png')\n"
            stext += "  for i in "+str(antNames)+":\n    aU.getTPSampling(vis = '"+msName+"', antenna = i, showplot = True, plotfile = '"+msName+".sampling.'+i+'.png')"

            sfsdr.addReducScriptStep(f1, mystepdict, "listobs", stext, mystepindent)

            stext = doAprioriFlagging(msName, flagAutoCorr=False, flagCalIntents=False, valueMaps=valueMaps)
            sfsdr.addReducScriptStep(f1, mystepdict, "A priori flagging", stext, mystepindent)

            if mycasaversion < '5.0':

                stext = "for i in "+str(antNames)+":\n  os.system('rm -Rf "+msName+".'+i+'*')\n\n"
                stext += "sdsave(infile = '"+msName+"',\n  splitant = True,\n  outfile = '"+msName+".asap',\n  overwrite = True)\n\n"
                sfsdr.addReducScriptStep(f1, mystepdict, "Split by antenna", stext, mystepindent)

                asapNames = [msName+'.'+i+'.asap' for i in antNames]
                asapNames = sorted(asapNames)

                stext = ''
                for i in asapNames:
                    stext += "os.system('rm -Rf "+i+".sdlist')\n"
                    stext += "sdlist(infile = '"+i+"',\n  outfile = '"+i+".sdlist')\n\n"
                sfsdr.addReducScriptStep(f1, mystepdict, "sdlist", stext, mystepindent)

            else:

                asapNames = [msName]

            tsysCalTableName = []
            stext = SDdoFillTsysSolutions(asapNames, msName=msName, doplot=True, tsysCalTableName=tsysCalTableName)

            sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the Tsys cal table", stext, mystepindent)

            skyCalTableName = []

            if mycasaversion >= '5.0':

                if step in ['SDcalibLine', 'SDscience']:
                    stext = SDdoFillTsysSolutions(asapNames, msName=msName, doplot=True, tsysCalTableName=skyCalTableName, sky=True, calmode='ps')
                else:
                    stext = SDdoFillTsysSolutions(asapNames, msName=msName, doplot=True, tsysCalTableName=skyCalTableName, sky=True, calmode='otfraster')
                sfsdr.addReducScriptStep(f1, mystepdict, "Generation of the Sky cal table", stext, mystepindent)

            spwInfo = sfsdr.getSpwInfo(msName, caching=True)
            if msName in valueMaps.keys():
                vm = valueMaps[msName]
                print('Using canned ValueMap.')
            else:
                vm = aU.ValueMapping(msName)
                valueMaps[msName] = vm

            jyCalTableName = []

            if step in ['SDcalibLine', 'SDscience'] and mycasaversion >= '6.4.3': # add K to Jy conversion step

                jyCalTableName.append(msName+'.jy')

                stext = "gencal(vis='"+msName+"',\n"
                stext += "   caltype='jyperk',\n"   
                stext += "   caltable='"+jyCalTableName[0]+"',\n"
                stext += "   spw = '"+','.join([str(j) for j in sorted(spwInfo.keys())])+"')\n"
                sfsdr.addReducScriptStep(f1, mystepdict, "Create caltable to convert the science target units from Kelvin to Jansky", stext, mystepindent)


            chanFlags = {}

            for i in sorted(spwInfo.keys()):
                if vm.spwInfo[i]['bandwidth'] > 1875000000:
                    chanEdge = int((vm.spwInfo[i]['numChannels'] - (1875000000. / vm.spwInfo[i]['bandwidth']) * vm.spwInfo[i]['numChannels']) / 2.)
                    maskflag = str([[0, chanEdge-1], [vm.spwInfo[i]['numChannels']-chanEdge, vm.spwInfo[i]['numChannels']-1]])
                    if maskflag not in list(chanFlags.keys()):
                        chanFlags[maskflag] = []
                    chanFlags[maskflag].append(i)

            stext = ''

            for i in asapNames:
                for maskflag in list(chanFlags.keys()):
                    maskflag1 = eval(maskflag)
                    maskflag2 = []
                    for j in range(len(maskflag1)):
                        maskflag2.append('~'.join([str(k) for k in maskflag1[j]]))
                    maskflag2 = ';'.join(maskflag2)
                    spwmaskflag = []
                    for j in chanFlags[maskflag]:
                        spwmaskflag.append(str(j)+':'+maskflag2)
                    spwmaskflag = ','.join(spwmaskflag)
                    if mycasaversion >= '5.0':
                        stext += "flagdata(vis = '"+i+"',\n  mode = 'manual',\n  spw = '"+spwmaskflag+"')\n\n"
                    else:
                        stext += "sdflag(infile = '"+i+"',\n  mode = 'manual',\n  spw = '"+spwmaskflag+"',\n  overwrite = True)\n\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Do initial flagging", stext, mystepindent)


            mytb.open(msName+'/OBSERVATION')
            obsTimeRange = mytb.getcol('TIME_RANGE')
            mytb.close()
            obsTimeStart = ((obsTimeRange[0]/86400.0)+2400000.5-2440587.5)*86400.0
            obsTimeStart = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTimeStart[0])))

            if step in ['SDcalibLine', 'SDscience']:

                if mycasaversion < '5.0':
                    stext = SDdoCalibration(asapNames, msName=msName, calmode='ps', tsysCalTableName=tsysCalTableName[0])
                    sfsdr.addReducScriptStep(f1, mystepdict, "Calibration of the data into Kelvins", stext, mystepindent)
                elif mycasaversion < '6.4.3':
                    stext = SDdoCalibration(asapNames, msName=msName, tsysCalTableName=tsysCalTableName[0], skyCalTableName=skyCalTableName[0])
                    sfsdr.addReducScriptStep(f1, mystepdict, "Calibration of the data into Kelvins", stext, mystepindent)
                else:
                    stext = SDdoCalibration(asapNames, msName=msName, tsysCalTableName=tsysCalTableName[0], skyCalTableName=skyCalTableName[0],
                                            jyCalTableName=jyCalTableName[0])
                    sfsdr.addReducScriptStep(f1, mystepdict, "Calibration of the data into Janskys", stext, mystepindent)

                    stext = SDdoAtmCor(msName=msName, jyCalTableName=jyCalTableName[0])
                    sfsdr.addReducScriptStep(f1, mystepdict, "Correction of residual atmospheric features (sdatmcor)", stext, mystepindent)


                if mycasaversion < '5.0':

                    asapNames = [i+'.cal' for i in asapNames]

                    if obsTimeStart < '2015-10-01T00:00:00':

                        stext = ''
                        for i in asapNames:
                            stext += "os.system('rm -Rf "+i+".nlc')\n\n"
                            stext += "sdscale(infile = '"+i+"',\n  outfile = '"+i+".nlc',\n  factor = 1.25)\n\n"
                        sfsdr.addReducScriptStep(f1, mystepdict, "Application of non-linearity correction factor", stext, mystepindent)

                        asapNames = [i+'.nlc' for i in asapNames]

                stext = SDdoBaselineSubtraction(asapNames, msName=msName)
                sfsdr.addReducScriptStep(f1, mystepdict, "Subtracting the baseline", stext, mystepindent)

                asapNames = [i+'.bl' for i in asapNames]

            else:

                if mycasaversion < '5.0':
                    stext = SDdoCalibration(asapNames, msName=msName, calmode='otfraster', tsysCalTableName=tsysCalTableName[0])
                else:
                    stext = SDdoCalibration(asapNames, msName=msName, tsysCalTableName=tsysCalTableName[0], skyCalTableName=skyCalTableName[0])
                sfsdr.addReducScriptStep(f1, mystepdict, "Calibration of the data into Kelvins", stext, mystepindent)

                if mycasaversion < '5.0':

                    asapNames = [i+'.cal' for i in asapNames]

                    if obsTimeStart < '2015-10-01T00:00:00':

                        stext = ''
                        for i in asapNames:
                            stext += "os.system('rm -Rf "+i+".nlc')\n\n"
                            stext += "sdscale(infile = '"+i+"',\n  outfile = '"+i+".nlc',\n  factor = 1.25)\n\n"
                        sfsdr.addReducScriptStep(f1, mystepdict, "Application of non-linearity correction factor", stext, mystepindent)

                        asapNames = [i+'.nlc' for i in asapNames]

            spwInfo = sfsdr.getSpwInfo(msName, caching=True)
            spwIds = sorted(spwInfo.keys())
            spwIds = [str(i) for i in spwIds]

            if mycasaversion < '5.0':

                stext = ''
                for i in asapNames:
                    stext += "os.system('rm -Rf "+i+".ms')\n\n"
                    if mycasaversion >= '4.2.2':
                        stext += "sdsave(infile = '"+i+"',\n  outfile = '"+i+".ms',\n  spw = '"+','.join(spwIds)+"',\n  outform = 'MS2')\n\n"
                    else:
                        stext += "sdsave(infile = '"+i+"',\n  outfile = '"+i+".ms',\n  outform = 'MS2')\n\n"
                sfsdr.addReducScriptStep(f1, mystepdict, "Converting ASAP -> MS", stext, mystepindent)

                asapNames = [i+'.ms' for i in asapNames]


            stext = ''

            if len(asapNames) > 1:
                stext += "os.system('rm -Rf "+msName+".cal')\n\n"
                stext += "concat(vis = [ \\\n    '"+"', \\\n    '".join(asapNames)+"' ], \\\n  concatvis = '"+msName+".cal')\n\n"
            else:
                stext += "os.system('rm -Rf "+msName+".cal')\n\n"
                stext += "os.system('cp -Rf "+asapNames[0]+" "+msName+".cal')\n\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Split and concatenation", stext, mystepindent)

            if step == 'SDscience' and mycasaversion < '6.4.3':

                jyperk = sfsdr.getJyPerK(msName, interactive=True)

                f2 = open('jyperk.txt', 'w')
                pprint.pprint(jyperk, stream=f2, indent=2)
                f2.close()

                f2 = open('jyperk.txt')
                jyperk = f2.read()
                f2.close()

                stext = "jyperk = \\\n" + jyperk + "\n"

                if mycasaversion >= '5.0':
                    blspwmap = {}
                    for i in spwIds: blspwmap[str(i)] = str(spwIds.index(i))
                    stext += "blspwmap = "+str(blspwmap)+"\n"

                stext += "os.system('rm -Rf "+msName+".cal.jy')\n"
                stext += "os.system('cp -Rf "+msName+".cal "+msName+".cal.jy')\n\n"

                stext += "for ant in jyperk.keys():\n"
                stext += "  for spw in jyperk[ant].keys():\n"

                if mycasaversion < '5.0':
                    stext += "    scaleAutocorr(vis='"+msName+".cal.jy', scale=jyperk[ant][spw]['mean'], antenna=ant, spw=spw)\n"
                else:
                    stext += "    scaleAutocorr(vis='"+msName+".cal.jy', scale=jyperk[ant][spw]['mean'], antenna=ant, spw=int(blspwmap[str(spw)]))\n"

                sfsdr.addReducScriptStep(f1, mystepdict, "Convert the Science Target Units from Kelvin to Jansky", stext, mystepindent)

            ###############################################

            if step in ['SDampcal']:

                print('Preparing imaging code ...')

                imagingParams[msName] = {}

                imagingParams[msName]['spwIds'] = spwIds
                print("running au.getTPSampling('%s', showplot=False, pickFirstRaster=True)" % (msName))
                xSampling, ySampling, maxsize = aU.getTPSampling(msName, showplot=False, pickFirstRaster=True)
                imagingParams[msName]['maxsize'] = float(maxsize)
                mymsmd = msmdtool()
                mymsmd.open(msName)

                for i in sorted(spwInfo.keys()):

                    imagingParams[msName][i] = {}

                    freq = mymsmd.meanfreq(i)
                    imagingParams[msName][i]['freq'] = freq

                    theorybeam = aU.primaryBeamArcsec(frequency=freq*1e-9, fwhmfactor=1.13, diameter=12)
                    imagingParams[msName][i]['theorybeam'] = theorybeam

                    minor, major, fwhmsfBeam, sfbeam = aU.sfBeam(frequency=freq*1e-9, pixelsize=theorybeam/9.0, convsupport=6, img=None, stokes='both', xSamplingArcsec=xSampling, ySamplingArcsec=ySampling, fwhmfactor=1.13, diameter=12)
                    imagingParams[msName][i]['sfbeam'] = sfbeam

                fieldId = mymsmd.fieldsforintent('OBSERVE_TARGET#ON_SOURCE')
                if len(fieldId) != 1: 
                    casalog.post('ERROR: UNEXPECTED NUMBER OF FIELDS.', 'SEVERE')
                    return False
                fieldId = fieldId[0]
                imagingParams[msName]['fieldId'] = fieldId

                fieldName = mymsmd.namesforfields(fieldId)
                if len(fieldName) != 1: 
                    casalog.post('ERROR: UNEXPECTED NUMBER OF FIELDS.', 'SEVERE')
                    return False
                fieldName = fieldName[0]
                imagingParams[msName]['fieldName'] = fieldName

                mymsmd.close()

            if step in ['SDscience']:

                print('Preparing imaging code (SDscience) ...')
                ignoreOFF = True # meaning we will prepare a POINTING table for getTPSampling which has no off-source positions
                myqa = qatool()

                imagingParams[msName] = {}

                imagingParams[msName]['spwIds'] = spwIds

                ephemFieldNames = aU.getEphemerisFields(msName)

                mymsmd = msmdtool()
                mymsmd.open(msName)

                fieldIds = list(mymsmd.fieldsforintent('OBSERVE_TARGET#ON_SOURCE'))

                isMultiField=False
                if len(fieldIds) > 1:
                    isMultiField=True
                    casalog.post('NOTE: there are '+str(len(fieldIds))+' science fields, i.e. more than one.', 'WARN')
                elif len(fieldIds) == 0:
                    casalog.post('ERROR: There are no science fields.', 'SEVERE')
                    return False

                imagingParams[msName]['fieldIds'] = fieldIds
                imagingParams[msName]['maxsize'] = {}
                fieldNames = {}

                for i in fieldIds:
                    fieldName = mymsmd.namesforfields(i)[0]
                    fieldNames[i] = fieldName

                    scansToUse = np.intersect1d(mymsmd.scansforfield(i), mymsmd.scansforintent('OBSERVE_TARGET#ON_SOURCE'))
                    if len(scansToUse)==0:
                        casalog.post('No scans observing field '+str(i), 'WARN')
                        imagingParams[msName]['maxsize'][fieldName] = 0.
                        continue
                    scansToUseStr = str(scansToUse[0])
                    if not (fieldName in ephemFieldNames): # for ephem fields, use only first scan, otherwise all relevant
                        for k in scansToUse[1:]:
                            scansToUseStr+=','+str(k)

                    if ignoreOFF:
                        if isMultiField:
                            mytimes = mymsmd.timesforfield(i)
                            myt = myqa.quantity(v=mytimes[0], unitname='s') 
                            myt2 = myqa.quantity(v=mytimes[-1], unitname='s') 
                            mytimerange = myqa.time(myt, form='ymd')[0]+'~'+myqa.time(myt2, form='ymd')[0]                        
                            print('Splitting out time range '+mytimerange+' to determine sampling for field '+fieldName)
                            tmp_msname = 'tmp_field'+str(i)+'_'+msName
                            os.system('rm -rf '+tmp_msname)
                            mstransform(vis=msName, timerange=mytimerange, spw=sorted(spwInfo.keys())[0], outputvis=tmp_msname, 
                                        nchan=1, datacolumn='data', antenna='0&&0', scan=scansToUseStr)  
                        else:
                            # we save the time for splitting and modify the POINTING table in place but keep a copy
                            tmp_msname = msName 
                            tmp_orig_pointing = msName+'/originalPOINTING'
                            os.system('rm -rf '+tmp_orig_pointing)
                            os.system('cp -R '+msName+'/POINTING '+tmp_orig_pointing)

                        # load from POINTING columns TIME, from Main TIME and STATE_ID, from STATE: OBS_MODE
                        mytb.open(tmp_msname+'/POINTING')
                        poiTime = mytb.getcol('TIME')
                        poiInterval = mytb.getcell('INTERVAL',0) 
                        mytb.close()
                        mytb.open(tmp_msname)
                        mainTime = mytb.getcol('TIME')
                        mainStateId = mytb.getcol('STATE_ID')
                        mytb.close()
                        mytb.open(tmp_msname+'/STATE')
                        stateObsMode = mytb.getcol('OBS_MODE')
                        mytb.close()
                        # Loop over OBS_MODE and get state_ids for ON_SOURCE
                        onsource_stateIds = []
                        for k in range(len(stateObsMode)):
                            if stateObsMode[k] == 'OBSERVE_TARGET#ON_SOURCE':
                                onsource_stateIds.append(k)
                        print('   onsource_stateIds ', onsource_stateIds)
                        print('   Determining on-source subscan times ...')
                        # loop over mainTime and mainStateId to get beginning and end of each on_source sub-scan: subscan_start[], subscan_end[]
                        subscan_start = []
                        subscan_end = []
                        timeSafetyMargin = 1.0 # seconds
                        if timeSafetyMargin < poiInterval:
                            timeSafetyMargin = poiInterval
                        curr_sId = -1 # the current onsource state ID
                        for k, mT in enumerate(mainTime):
                            if curr_sId > 0:
                                if mainStateId[k] == curr_sId:
                                    continue
                                else: # we have reached the end of an onsource subscan
                                    subscan_end.append(mT-timeSafetyMargin)
                                    curr_sId = -1
                            else: # search for next 
                                for sId in onsource_stateIds:
                                    if mainStateId[k] == sId: # we are in an onsource subscan
                                        subscan_start.append(mT+timeSafetyMargin)
                                        curr_sId = sId
                                        break
                        if len(subscan_start)>len(subscan_end): 
                            # the last onsource subscan ended at the end of mainTime
                            subscan_end.append(mainTime[-1]-timeSafetyMargin)

                        # loop over pointing time and compile list of rows which are not in an on_source subscan
                        print('   Removing other subscan times from POINTING ...')
                        pRowsToBeDel = []
                        kstart = 0
                        for j, pT in enumerate(poiTime):
                            notFound = True
                            for k in range(kstart, len(subscan_start)):
                                if subscan_start[k]<=pT and pT<subscan_end[k]:
                                    notFound = False
                                    kstart = k # earlier subscans can now be ignored
                                    break
                            if notFound:
                                #print(j)
                                pRowsToBeDel.append(j)
                        # delete the list of rows from POINTING
                        mytb.open(tmp_msname+'/POINTING', nomodify=False)
                        print('   Deleting '+str(len(pRowsToBeDel))+' rows from POINTING table of '+tmp_msname)
                        mytb.removerows(pRowsToBeDel)
                        mytb.close()
                        samplingPlotName = msName+'.'+'sampling_field'+str(i)+'.png'

                        print("   running au.getTPSampling('"+tmp_msname+"', showplot=True, plotfile='"+samplingPlotName+"',pickFirstRaster=False, field='"+fieldName+"', scan='"+scansToUseStr+"')")
                        xSampling, ySampling, maxsize = aU.getTPSampling(tmp_msname, showplot=True, plotfile=samplingPlotName, pickFirstRaster=False, field=fieldName, scan=scansToUseStr)

                        if isMultiField:
                            os.system('rm -rf '+tmp_msname)
                        else:
                            print('Restoring POINTING table in '+msName)
                            os.system('rm -rf '+msName+'/POINTING')
                            os.system('mv '+tmp_orig_pointing+' '+msName+'/POINTING')

                    else: # don't ignore off
                        print("running au.getTPSampling("+msName+", showplot=False, pickFirstRaster=False, field='"+fieldName+"', scan="+scansToUseStr+")")
                        xSampling, ySampling, maxsize = aU.getTPSampling(msName, showplot=False, pickFirstRaster=False, field=fieldName, scan=scansToUseStr)

                    print("   ... found maxsize (arcsec) for field "+str(i)+" to be ", float(maxsize))
                    imagingParams[msName]['maxsize'][fieldName] = float(maxsize)

                for i in sorted(spwInfo.keys()):
                    imagingParams[msName][i] = {}

                    freq = mymsmd.meanfreq(i)
                    imagingParams[msName][i]['freq'] = freq

                    theorybeam = aU.primaryBeamArcsec(frequency=freq*1e-9, fwhmfactor=1.13, diameter=12)
                    imagingParams[msName][i]['theorybeam'] = theorybeam

                    minor, major, fwhmsfBeam, sfbeam = aU.sfBeam(frequency=freq*1e-9, pixelsize=theorybeam/9.0, convsupport=6, img=None, stokes='both', xSamplingArcsec=xSampling, ySamplingArcsec=ySampling, fwhmfactor=1.13, diameter=12)
                    imagingParams[msName][i]['sfbeam'] = sfbeam

                mymsmd.close()


            if step == 'SDampcal':

                stext = "# the values below were calculated assuming fwhmfactor = 1.13\n\n"
                stext += "maxsize = "+str(imagingParams[msName]['maxsize'])+"\n\n"

                stext += "theorybeam = {}\n"
                for i in sorted(spwInfo.keys()):
                    stext += "theorybeam['"+str(i)+"'] = "+str(imagingParams[msName][i]['theorybeam'])+" # mean freq = "+str(imagingParams[msName][i]['freq']*1e-9)+"\n"

                setjyModels = ['Venus', 'Mars', 'Jupiter', 'Uranus', 'Neptune', 'Pluto', 'Io', 'Europa', 'Ganymede', 'Callisto', 'Titan', 'Triton', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Victoria', 'Davida']

                doPlanet = 0
                for j in range(len(setjyModels)):
                    if re.search(setjyModels[j], fieldName, re.IGNORECASE) is not None:
                        doPlanet = 1
                        break

                if doPlanet == 1:
                    supportedSSOTPampCals = ['mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
                    if fieldName.lower() not in supportedSSOTPampCals: 
                        casalog.post('ERROR: SSO NOT SUPPORTED AS TP AMP CAL.', 'SEVERE')

                stext += "\nfor spw in "+str(spwIds)+":\n\n"
                stext += "  cell = theorybeam[spw]/9.0\n"
                stext += "  imsize = int(1.5*maxsize/cell)\n"
                stext += "  imsize += (imsize % 2)\n\n"
                stext += "  for ant in "+str(antNames)+":\n\n"
                stext += "    sdimaging(infiles = '"+msName+".cal',\n"
                stext += "      field = '"+fieldName+"',\n"
                stext += "      spw = spw,\n"
                stext += "      antenna = ant,\n"
                stext += "      nchan = 1,\n"
                stext += "      mode = 'channel',\n"
                stext += "      width = '4080',\n"
                stext += "      gridfunction = 'SF',\n"
                stext += "      convsupport = 6,\n"

                if doPlanet == 1:
                    stext += "      ephemsrcname = '"+fieldName+"',\n"
                else:
                    stext += "      phasecenter = "+str(fieldId)+",\n"

                stext += "      imsize = imsize,\n"
                stext += "      cell = str(cell)+'arcsec',\n"
                stext += "      overwrite = True,\n"
                stext += "      outfile = '"+msName+".cal.%s.spw%s.image' % (ant, spw))\n"

                sfsdr.addReducScriptStep(f1, mystepdict, "Imaging", stext, mystepindent)

                stext = "srcflux = {}\n"

                for i in sorted(spwInfo.keys()):

                    if doPlanet == 1:
                        planetInfo = aU.planetFlux(vis=msName, spw=i)
                        srcflux = planetInfo['fluxDensity']
                        srcsize = (planetInfo['majorAxis']*planetInfo['minorAxis'])**0.5
                        spwfreq = planetInfo['meanFrequency']
                    else:
                        qsoInfo = aU.getALMAFluxForMS(msName, field=fieldName, spw=str(i), useCalibratorService=useCalibratorService,
                                                      calibratorServiceURL=calibratorServiceURL)
                        srcflux = qsoInfo[fieldName]['fluxDensity']
                        spwfreq = qsoInfo[fieldName]['frequency']

                    stext += "srcflux['"+str(i)+"'] = "+str(srcflux)+" # mean freq = "+str(spwfreq*1e-9)+"\n"

                if doPlanet == 1:
                    stext += "\nsrcsize = "+str(srcsize)+"\n"

                stext += "\njyperk = {}\n"

                stext += "\nfor ant in "+str(antNames)+":\n\n"
                stext += "  jyperk[ant] = {}\n\n"
                stext += "  for spw in "+str(spwIds)+":\n\n"
                stext += "    if os.path.exists('"+msName+".cal.%s.spw%s.image' % (ant, spw)):\n\n"
                stext += "      peak = imstat('"+msName+".cal.%s.spw%s.image' % (ant, spw))['max'][0]\n\n"

                if doPlanet == 1:

                    stext += "      fwhm = aU.getfwhm2('"+msName+".cal.%s.spw%s.image' % (ant, spw))\n"
                    stext += "      print('Apparent FWHM (inc. gridding convolution) is %.2f arcsec' % fwhm)\n\n"
                    stext += "      deconvfwhm = aU.deconvolveDiskFromBeam(fwhm, srcsize)\n"
                    stext += "      print('FWHM deconvolved for planet size is %.2f arcsec' % deconvfwhm)\n\n"
                    stext += "      # correction factor for dilution due to planet size\n"
                    stext += "      srcsizecorr = (deconvfwhm/fwhm)**2\n\n"
                    stext += "      jyperk[ant][spw] = srcflux[spw] / peak * srcsizecorr"

                else:

                    stext += "      jyperk[ant][spw] = srcflux[spw] / peak"

                sfsdr.addReducScriptStep(f1, mystepdict, "Determination of the Jy/K factors", stext, mystepindent)

                stext = "asdm = '"+re.findall('uid___[a-zA-Z0-9]+_[a-zA-Z0-9]+_[a-zA-Z0-9]+', msName)[0]+"'\n"

                mytb.open(msName+'/OBSERVATION')
                obsTimeRange = mytb.getcol('TIME_RANGE')
                obsTime = (obsTimeRange[0]+obsTimeRange[1])/2.0
                obsTime = ((obsTime/86400.0)+2400000.5-2440587.5)*86400.0
                obsTime = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTime[0])))
                mytb.close()

                stext += "date = '"+obsTime+"'\n"

                stext += "ampcal = '"+fieldName+"'\n"
                mymsmd = msmdtool()
                mymsmd.open(msName)

                spwbb = {}
                spwfreq = {}
                spwbw = {}

                for i in sorted(spwInfo.keys()):
                    spwbb[str(i)] = str(mymsmd.baseband(int(i)))
                    spwfreq[str(i)] = str(mymsmd.meanfreq(int(i)))
                    spwbw[str(i)] = str(mymsmd.bandwidths(int(i)))

                stext += "band = '"+str(aU.getBand(spwfreq[list(spwfreq.keys())[0]]))+"'\n"
                stext += "bb = "+str(spwbb)+" # spw baseband number\n"
                stext += "freq = "+str(spwfreq)+" # spw mean frequency\n"
                stext += "bw = "+str(spwbw)+" # spw bandwidth\n"

                scanList = mymsmd.scansforintent('OBSERVE_TARGET#ON_SOURCE')

                mymsmd.close()

                weatherInfo = aU.getWeather(msName, scan = scanList.tolist(), getSolarDirection=False)

                stext += "elev = '"+str(weatherInfo[0]['elevation'])+"' # mean elevation\n"
                stext += "temp = '"+str(weatherInfo[0]['temperature'])+"' # mean temperature\n\n"

                stext += "f = open('"+msName+".cal.jyperk.txt', 'w')\n\n"

                stext += "for ant in "+str(antNames)+":\n"
                stext += "  for spw in "+str(spwIds)+":\n\n"
                stext += "    if ant in jyperk.keys():\n"
                stext += "      if spw in jyperk[ant].keys():\n"
                stext += "        f.write('\\t'.join([asdm, ant, spw, str(jyperk[ant][spw]), date, ampcal, band, bb[spw], freq[spw], bw[spw], elev, temp])+'\\n')\n"

                stext += "\nf.close()\n"

                sfsdr.addReducScriptStep(f1, mystepdict, "Writing out the Jy/K factors", stext, mystepindent)

            #endif step = 'SDampcal'

            sfsdr.prependReducScriptHeader(f1, mystepdict, "Created using "+version(), mystepindent)
            print("Wrote ", f1name)
            f1.close()

        ###################################
        
        if step == 'SDscience':

            os.chdir(currDir)

            mystepdict = {}
            mystepindent = "  "
            f1name = 'scriptForSDimaging.py'
            f1 = open(f1name, 'w')

            if mycasaversion < '5.1':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', casaVersion) == None:", file=f1)
            elif mycasaversion < '5.9':
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in cu.version().tolist()[:-1]])) == None:", file=f1)
            else:
                print("if re.search('^"+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"', '.'.join([str(i) for i in casalith.version()[:-1]])) == None:", file=f1)
            print("  sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: "+re.findall('^[0-9]+.[0-9]+.[0-9]+', mycasaversion)[0]+"')\n\n", file=f1)

            stext = 'msNames = [ \\\n'
            for i in range(len(msNames)):
                if mycasaversion < '6.4.4':
                    stext += "'"+msNames[i]+".cal.jy', \\\n"
                else:
                    stext += "'"+msNames[i]+".cal', \\\n"
            stext += ']\n\n'

            tos = []
            stext += "# Time on source\n"
            for i in range(len(msNames)):
                tos.append(aU.timeOnSource(msNames[i])['minutes_on_science'])
                stext += "# "+msNames[i]+" = "+str(round(tos[i], 1))+" min\n"
            stext += "# Total = "+str(round(np.sum(tos), 1))+" min"

            sfsdr.addReducScriptStep(f1, mystepdict, "Define the calibrated datasets", stext, mystepindent)

            ############################################

            splitMS = False
            for i in range(len(msNames)):
                if sorted(imagingParams[msNames[i]]['spwIds']) != sorted(imagingParams[msNames[0]]['spwIds']):
                    print('WARNING: the science spw indices are not the same across all EBs, I will need to run split.')
                    splitMS = True
                    break

            if splitMS == True:

                stext = ''
                for i in range(len(msNames)):

                    stext += "os.system('rm -Rf "+msNames[i]+".cal.jy.split')\n\n"
                    stext += "split(vis = '"+msNames[i]+".cal.jy',\n"
                    stext += "  outputvis = '"+msNames[i]+".cal.jy.split',\n"
                    stext += "  datacolumn = 'all',\n"
                    stext += "  spw = '"+','.join(sorted(imagingParams[msNames[i]]['spwIds']))+"')\n\n"

                    spwIds1 = sorted(imagingParams[msNames[i]]['spwIds'])
                    spwIds1 = [int(j) for j in spwIds1]
                    for j in range(len(spwIds1)):
                        imagingParams[msNames[i]][j] = imagingParams[msNames[i]][spwIds1[j]]
                        imagingParams[msNames[i]].pop(spwIds1[j])

                    spwIds1 = list(range(len(imagingParams[msNames[i]]['spwIds'])))
                    spwIds1 = [str(j) for j in spwIds1]
                    imagingParams[msNames[i]]['spwIds'] = spwIds1

                stext += 'msNames = [ \\\n'
                for i in range(len(msNames)):
                    stext += "'"+msNames[i]+".cal.jy.split', \\\n"
                stext += ']\n\n'

                sfsdr.addReducScriptStep(f1, mystepdict, "Split the science spectral windows", stext, mystepindent)

            #######################################################


            myspwids = sorted(imagingParams[msNames[0]]['spwIds'])

            stext = "# only the spectral windows listed below will be imaged\n"

            stext += "spwIds = "+str(myspwids)+"\n\n"

            if mycasaversion >= '5.0':

                blspwmap = {}
                for i in myspwids: blspwmap[i] = str(myspwids.index(i))
                stext += "blspwmap = "+str(blspwmap)+"\n"

            # image name endings 
            myimagename = {}            
            for i in myspwids:
                myimagename[i] = '_sci.spw'+str(i)+'.cube.I.manual'

            stext += "\nmyimagenames = {"
            for i in myspwids:
                stext += "'"+str(i)+"': '"+myimagename[i]+"',\n                "
            stext += "}\n"

            # science fields
            stext += "\nmyfields = {"
            for i in fieldNames.keys():
                stext += "'"+str(i)+"': '"+fieldNames[i]+"',\n            "
            stext += "}\n"

            # maxsizes for the science fields
            stext += "\nmymaxsizes = {"
            for i in fieldNames.keys():
                stext += "'"+str(i)+"': "+str(imagingParams[msNames[0]]['maxsize'][fieldNames[i]])+",\n              "
            stext += "}\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Define the imaging parameters", stext, mystepindent)

            #######################################################

            haveEphemTarget = False
            for i in fieldNames.keys():
                if fieldNames[i] in ephemFieldNames:
                    haveEphemTarget = True
                    break

            if mycasaversion < '6.6' and not haveEphemTarget:
                stext = ''
                sdimgPrefix = '' # use oldfashioned sdimaging
            else:
                stext = "import os\n\n"
                sdimgPrefix = 't' # use tsdimaging
                
            stext += "# the values below were calculated assuming fwhmfactor = 1.13\n\n"

            stext += "theorybeam = {}\n"

            for i in myspwids:
                stext += "theorybeam['"+str(i)+"'] = "+str(imagingParams[msNames[0]][int(i)]['theorybeam'])+" # mean freq = "+str(imagingParams[msNames[0]][int(i)]['freq']*1e-9)+"\n"

            stext += "\n# the values below were calculated assuming cell = theorybeam[spw]/9.0\n"
            stext += "sfbeam = {}\n"

            for i in myspwids:
                stext += "sfbeam['"+str(i)+"'] = "+str(imagingParams[msNames[0]][int(i)]['sfbeam'])+" # mean freq = "+str(imagingParams[msNames[0]][int(i)]['freq']*1e-9)+"\n"

            # improvement to help with cases like SACM-576 
            limit_numebs = 9
            if len(msNames) > limit_numebs:
                stext += "\n# This dataset includes more than "+str(limit_numebs)+" EBs which may cause memory problems for sdimaging\n"
                stext += "# when using a list of MSs as parameter 'infiles'. Working with concatenated dataset instead.\n"
                
                stext += "\nconcat(vis = msNames, concatvis = 'concat_"+str(len(msNames))+"EBs.ms')\n"

            stext += "\nfor myfieldid in myfields.keys():\n\n"

            stext += "  maxsize = mymaxsizes[myfieldid]\n"

            stext += "\n  for spw in spwIds:\n\n"

            stext += "    cell = theorybeam[spw]/9.0\n"
            stext += "    imsize = int(1.5*maxsize/cell)\n\n"
            stext += "    imsize += (imsize % 2)\n\n"


            if len(msNames) > limit_numebs:
                stext += "    "+sdimgPrefix+"sdimaging(infiles = 'concat_"+str(len(msNames))+"EBs.ms',\n"
            else:
                stext += "    "+sdimgPrefix+"sdimaging(infiles = msNames,\n"


            stext += "      field = myfields[myfieldid],\n"

            if mycasaversion < '5.0':
                stext += "      spw = spw,\n"
            else:
                stext += "      spw = blspwmap[spw],\n"

            stext += "      mode = 'channel',\n"
            stext += "      outframe = 'lsrk',\n"   # requires lower-case prior to fix for CAS-12820
            stext += "      gridfunction = 'SF',\n"
            stext += "      convsupport = 6,\n"
            if haveEphemTarget:
                stext += "      phasecenter = 'TRACKFIELD',\n"
                stext += "      specmode = 'cubesource',\n"
            else:
                stext += "      phasecenter = int(myfieldid),\n"

            stext += "      imsize = imsize,\n"
            stext += "      cell = str(cell)+'arcsec',\n"
            stext += "      overwrite = True,\n"
            stext += "      outfile = myfields[myfieldid]+myimagenames[spw])\n\n"

            if sdimgPrefix == 't': # make up for tsdiaging file naming issue
                stext += "    os.system('mv '+myfields[myfieldid]+myimagenames[spw]+'.image '+myfields[myfieldid]+myimagenames[spw])\n\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Image the Science Target", stext, mystepindent)

            stext = "\nfor myfieldid in myfields.keys():\n"
            stext += "\n  for spw in spwIds:\n\n"
            stext += "    imhead(imagename = myfields[myfieldid]+myimagenames[spw],\n"
            stext += "      mode = 'put',\n"
            stext += "      hdkey = 'bunit',\n"
            stext += "      hdvalue = 'Jy/beam')\n\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Correct the brightness unit in the image header", stext, mystepindent)

            stext = "\nfor myfieldid in myfields.keys():\n"
            stext += "\n  for spw in spwIds:\n\n"
            stext += "    myia = iatool()\n"
            stext += "    myia.open(myfields[myfieldid]+myimagenames[spw])\n"
            stext += "    myia.setrestoringbeam(major = str(sfbeam[spw])+'arcsec', minor = str(sfbeam[spw])+'arcsec', pa = '0deg')\n"
            stext += "    myia.done()\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Add Restoring Beam Header Information to the Science Image", stext, mystepindent)

            stext = "\nfor myfieldid in myfields.keys():\n"
            stext += "\n  for spw in spwIds:\n\n"
            stext += "    exportfits(imagename = myfields[myfieldid]+myimagenames[spw],\n"
            stext += "      fitsimage = myfields[myfieldid]+myimagenames[spw]+'.fits')\n\n"

            sfsdr.addReducScriptStep(f1, mystepdict, "Export images to fits", stext, mystepindent)

            sfsdr.prependReducScriptHeader(f1, mystepdict, "SD Imaging\n# Created using "+version(), mystepindent)

            print("Wrote ", f1name)
            f1.close()

    return True

    # end of generateReducScript()

############################

def doAprioriFlagging(msName, flagAutoCorr=True, flagCalIntents=True, valueMaps={}):
    """Generate code for the Apriori Flagging step of a calibration script."""

    print('\n*** doAprioriFlagging ***')

    casaCmd = ''

    print('Gathering information ...')

    if flagCalIntents == True:

        intentsToFlag = ['POINTING', 'FOCUS', 'SIDEBAND_RATIO', 'ATMOSPHERE']

        if msName in valueMaps.keys():
            vm = valueMaps[msName]
            print('Using canned ValueMap.')
        else:
            vm = aU.ValueMapping(msName)
            valueMaps[msName] = vm

        fullIntentList = vm.uniqueIntents

        scanIntentList = []
        for i in intentsToFlag:
            for j in fullIntentList:
                if re.search(i, j) is not None:
                    scanIntentList.append('*'+i+'*')
                    break

        scanIntentList = ','.join(scanIntentList)

    mytb = aU.createCasaTool(tbtool)
    if flagAutoCorr == True:

        mytb.open(msName+'/DATA_DESCRIPTION')
        spwIds = mytb.getcol('SPECTRAL_WINDOW_ID')
        mytb.close()

        mytb.open(msName+'/PROCESSOR')
        procType = mytb.getcol('TYPE')
        mytb.close()

        procType1 = np.where(procType == 'RADIOMETER')[0]

        spwIds1 = []

        mytb.open(msName)

        for i in procType1:
            tb1 = mytb.query('PROCESSOR_ID == '+str(i))
            dataDescIds1 = tb1.getcol('DATA_DESC_ID')
            dataDescIds1 = np.unique(dataDescIds1)
            for j in dataDescIds1:
                spwIds1.append(spwIds[j])
            tb1.close() 

        mytb.close()

        spwIds1 = [i for i in range(len(spwIds)) if i not in spwIds1]

        if len(spwIds1) > 1:
            j0 = 0
            spwIds2 = str(spwIds1[j0])
            for j in range(len(spwIds1)-1):
                if spwIds1[j+1] == spwIds1[j]+1: continue
                spwIds2 = spwIds2 + '~' + str(spwIds1[j])
                j0 = j+1
                spwIds2 = spwIds2 + ',' + str(spwIds1[j0])
            spwIds2 = spwIds2 + '~' + str(spwIds1[j+1])
        else:
            spwIds2 = str(spwIds1[0])

    if flagAutoCorr == True:
        if re.search('^3.3', aU.getCasaVersion()) is not None:
            casaCmd = casaCmd + "flagdata(vis = '"+msName+"',\n"
            casaCmd = casaCmd + "  mode = 'manualflag',\n"
            casaCmd = casaCmd + "  spw = '"+spwIds2+"',\n"
            casaCmd = casaCmd + "  autocorr = True,\n"
            casaCmd = casaCmd + "  flagbackup = False)\n\n"
        elif aU.getCasaVersion() >= '4.1.0':
            casaCmd = casaCmd + "flagdata(vis = '"+msName+"',\n"
            casaCmd = casaCmd + "  mode = 'manual',\n"
            casaCmd = casaCmd + "  spw = '"+spwIds2+"',\n"
            casaCmd = casaCmd + "  autocorr = True,\n"
            casaCmd = casaCmd + "  flagbackup = False)\n\n"
        else:
            casaCmd = casaCmd + "tflagdata(vis = '"+msName+"',\n"
            casaCmd = casaCmd + "  mode = 'manual',\n"
            casaCmd = casaCmd + "  spw = '"+spwIds2+"',\n"
            casaCmd = casaCmd + "  autocorr = True,\n"
            casaCmd = casaCmd + "  flagbackup = False)\n\n"

    if flagCalIntents == True:

        if re.search('^3.3', aU.getCasaVersion()) is not None:
            casaCmd = casaCmd + "flagdata(vis = '"+msName+"',\n"
            casaCmd = casaCmd + "  mode = 'manualflag',\n"
            casaCmd = casaCmd + "  intent = '"+scanIntentList+"',\n"
            casaCmd = casaCmd + "  flagbackup = False)\n"
        elif aU.getCasaVersion() >= '4.1.0':
            casaCmd = casaCmd + "flagdata(vis = '"+msName+"',\n"
            casaCmd = casaCmd + "  mode = 'manual',\n"
            casaCmd = casaCmd + "  intent = '"+scanIntentList+"',\n"
            casaCmd = casaCmd + "  flagbackup = False)\n"
        else:
            casaCmd = casaCmd + "tflagdata(vis = '"+msName+"',\n"
            casaCmd = casaCmd + "  mode = 'manual',\n"
            casaCmd = casaCmd + "  intent = '"+scanIntentList+"',\n"
            casaCmd = casaCmd + "  flagbackup = False)\n"

    mytb.open(msName)
    tableNames = mytb.keywordnames()
    mytb.close()

    if 'FLAG_CMD' in tableNames:

        mytb.open(msName+'/FLAG_CMD')
        nFlagRows = mytb.nrows()
        mytb.close()

        if nFlagRows != 0:

            if re.search('^3.3', aU.getCasaVersion()) == None:
                casaCmd = casaCmd + "\nflagcmd(vis = '"+msName+"',\n"
                casaCmd = casaCmd + "  inpmode = 'table',\n"
                casaCmd = casaCmd + "  useapplied = True,\n"
                casaCmd = casaCmd + "  action = 'plot',\n"
                casaCmd = casaCmd + "  plotfile = '"+msName+".flagcmd.png')\n\n"
                casaCmd = casaCmd + "flagcmd(vis = '"+msName+"',\n"
                casaCmd = casaCmd + "  inpmode = 'table',\n"
                casaCmd = casaCmd + "  useapplied = True,\n"
                casaCmd = casaCmd + "  action = 'apply')\n"
            else:
                casaCmd = casaCmd + "\nflagcmd(vis = '"+msName+"',\n"
                casaCmd = casaCmd + "  flagmode = 'table',\n"
                casaCmd = casaCmd + "  optype = 'plot')\n\n"
                casaCmd = casaCmd + "flagcmd(vis = '"+msName+"',\n"
                casaCmd = casaCmd + "  flagmode = 'table',\n"
                casaCmd = casaCmd + "  optype = 'apply')\n"

    return casaCmd

####################

def doGenerateWVRCalTable(msName, calTableName=[], refant='', smooth=True, doplot=True, remcloud=False, valueMaps={}, isB2B=False, isBWSW=False):
    """Generate code for the wvrgcal step of a calibration script."""

    print('\n*** doGenerateWVRCalTable ***')

    print('Gathering information ...')

    if remcloud:
        if aU.getCasaSubversionRevision() < '35187':
            casalog.post('ERROR: remcloud option is only supported for CASA >= r35187','SEVERE')
            return False
        print('Will include call to remove_cloud ...')

    casaCmd = ''
    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/OBSERVATION')
    obsTimeRange = mytb.getcol('TIME_RANGE')
    mytb.close()
    obsTimeStart = ((obsTimeRange[0]/86400.0)+2400000.5-2440587.5)*86400.0
    obsTimeStart = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTimeStart[0])))
    if obsTimeStart > '2013-01-21T00:00:00':
        wvrTimeOffset = 0
    else:
        wvrTimeOffset = -1

    mytb.open(msName+'/ANTENNA')
    antNames = mytb.getcol('NAME')
    mytb.close()

    found = 0
    for i in antNames:
        if re.search('^[CP]M[0-9]+', i) == None: found = 1
    if found == 0:
        calTableName.append('')
        return casaCmd

    intentSources = sfsdr.getIntentsAndSourceNames(msName)

    sciSourceId = intentSources['OBSERVE_TARGET']['sourceid']
    sciSourceId1 = list(dict.fromkeys(sciSourceId).keys())
    if sciSourceId1[0] != '':
        sciSourceName = intentSources['OBSERVE_TARGET']['name'][sciSourceId.index(min(sciSourceId1))]
        sciSourceId = min(sciSourceId1)
        phaseCal = sfsdr.getPhaseCal(msName, valueMaps=valueMaps)

    if remcloud:
        if aU.getCasaVersion() < '5.9.9':
            casaCmd = casaCmd + "import recipes.remove_cloud as rc\n\n"
        else:
            casaCmd = casaCmd + "import casarecipes.remove_cloud as rc\n\n"

    casaCmd = casaCmd + "os.system('rm -rf %s.wvr') \n\n"%(msName)

    if remcloud:
        casaCmd = casaCmd + "os.system('rm -rf %s.cloud_offsets') \n\n"%(msName)

    casaCmd = casaCmd + "os.system('rm -rf %s.wvrgcal') \n\n"%(msName)

    sciSpwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET', caching=True)
    spwInfo = sciSpwInfo.copy()
    if isB2B or isBWSW:
        refSpwInfo = sfsdr.getSpwInfo(msName, intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
        spwInfo.update(refSpwInfo)

    integTime = []
    for i in sciSpwInfo: 
        integTime.append(sciSpwInfo[i]['integTime'])
    integTime = list(dict.fromkeys(integTime).keys())
    if len(integTime) != 1: 
        casaCmd = casaCmd + "# Warning: more than one integration time found on science data, I'm picking the lowest value. Please check this is right.\n\n"
    integTime = min(integTime)

    if remcloud:
        casaCmd = casaCmd + "rc.remove_cloud(vis='"+msName+"', offsetstable='"+msName+".cloud_offsets')\n\n"

    casaCmd = casaCmd + "mylogfile = casalog.logfile()\n"
    casaCmd = casaCmd + "casalog.setlogfile('"+msName+".wvrgcal')\n\n"
    casaCmd = casaCmd + "wvrgcal(vis = '"+msName+"',\n"

    if remcloud:
        casaCmd = casaCmd + "  offsetstable = '"+msName+".cloud_offsets',\n"

    casaCmd = casaCmd + "  caltable = '"+msName+".wvr',\n"

    spwIds = sorted(spwInfo.keys())
    casaCmd = casaCmd + "  spw = "+str(spwIds)+",\n"

    if smooth == True and integTime > 1.152:
        casaCmd = casaCmd + "  smooth = '"+str(integTime)+"s',\n"

    casaCmd = casaCmd + "  toffset = "+str(wvrTimeOffset)
    if sciSourceId1[0] != '':
        casaCmd = casaCmd + ",\n  tie = "+str([','.join([i, phaseCal[i]['phaseCalName']]) for i in phaseCal])+",\n"
        casaCmd = casaCmd + "  statsource = '"+sciSourceName+"')\n\n"
    else:
        casaCmd = casaCmd + ")\n\n"
    casaCmd = casaCmd + "casalog.setlogfile(mylogfile)\n\n"

    calTableName1 = msName+'.wvr'


    calTableName.append(calTableName1)

    if doplot and isB2B:
        sciSpw0 = sorted(sciSpwInfo.keys())[0] 
        refSpw0 = sorted(refSpwInfo.keys())[0] 
        casaCmd = casaCmd + "if applyonly != True: aU.plotWVRSolutions(caltable='%s', spw='%s', antenna='%s',\n" %(calTableName1,sciSpw0,refant)
        casaCmd = casaCmd + "  yrange=[-199,199],subplot=22, interactive=False,\n"
        casaCmd = casaCmd + "  figfile='%s') \n\n" %(calTableName1+'.plots/'+calTableName1.split('/')[-1]+'SignalSPWs')
        casaCmd = casaCmd + "if applyonly != True: aU.plotWVRSolutions(caltable='%s', spw='%s', antenna='%s',\n" %(calTableName1,refSpw0,refant)
        casaCmd = casaCmd + "  yrange=[-199,199],subplot=22, interactive=False,\n"
        casaCmd = casaCmd + "  figfile='%s') \n\n" %(calTableName1+'.plots/'+calTableName1.split('/')[-1]+'ReferenceSPWs')
        casaCmd = casaCmd + "#Note: If you see wraps in these plots, try changing yrange or unwrap=True \n"
        casaCmd = casaCmd + "#Note: If all plots look strange, it may be a bad WVR on the reference antenna.\n"
        casaCmd = casaCmd + "#      To check, you can set antenna='' to show all baselines.\n"

    elif doplot:
        sciSpw0 = sorted(sciSpwInfo.keys())[0] 
        casaCmd = casaCmd + "if applyonly != True: aU.plotWVRSolutions(caltable='%s', spw='%s', antenna='%s',\n" %(calTableName1,sciSpw0,refant)
        casaCmd = casaCmd + "  yrange=[-199,199],subplot=22, interactive=False,\n"
        casaCmd = casaCmd + "  figfile='%s') \n\n" %(calTableName1+'.plots/'+calTableName1.split('/')[-1])
        casaCmd = casaCmd + "#Note: If you see wraps in these plots, try changing yrange or unwrap=True \n"
        casaCmd = casaCmd + "#Note: If all plots look strange, it may be a bad WVR on the reference antenna.\n"
        casaCmd = casaCmd + "#      To check, you can set antenna='' to show all baselines.\n"



    return casaCmd

###################################

def doRunFixPlanets(msName):
    """Generate code for running fixplanets on fields with (0,0) coordinates"""

    print('\n*** doRunFixPlanets ***')

    print('Gathering information ...')

    fieldIds = sfsdr.getFieldsForFixPlanets(msName)

    if len(fieldIds) != 0:

        casaCmd = ''
        mytb = aU.createCasaTool(tbtool)
        mytb.open(msName+'/FIELD')
        fieldNames = mytb.getcol('NAME')
        mytb.close()

        fieldNames = ['%s' %fieldNames[i] for i in fieldIds]
        fieldNames = ','.join(fieldNames)
        fieldIds = ['%s' %i for i in fieldIds]
        fieldIds = ','.join(fieldIds)

        casaCmd = casaCmd + "fixplanets(vis = '"+msName+"',\n"
        casaCmd = casaCmd + "  field = '"+fieldIds+"', # "+fieldNames+"\n"
        casaCmd = casaCmd + "  fixuvw = True)\n"

        return casaCmd

####################################

def doGenerateTsysCalTable(msName, calTableName=[], doplot=True, isB2B=False):
    """Generate code for the Tsys table generation step of a calibration script."""

    print('\n*** doGenerateTsysCalTable ***')
    print('Gathering information ...')

    casaCmd = ''
    
    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/OBSERVATION')
    obsTimeRange = mytb.getcol('TIME_RANGE')
    mytb.close()
    obsTimeStart = ((obsTimeRange[0]/86400.0)+2400000.5-2440587.5)*86400.0
    obsTimeStart = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTimeStart[0])))

    sciSpwInfo = sfsdr.getSpwInfo(msName, caching=True)

    tsysNumChans = []
    tsysSpwInfo = sfsdr.getSpwInfo(msName, intent='CALIBRATE_ATMOSPHERE', caching=True)
    for i in tsysSpwInfo: 
        tsysNumChans.append(tsysSpwInfo[i]['numChans'])
    tsysNumChans = sorted(dict.fromkeys(tsysNumChans).keys())

    tsysNumChans = tsysNumChans[0]

    casaCmd = casaCmd + "os.system('rm -rf %s.tsys') \n"%(msName)  # Added by CLB

    casaCmd = casaCmd + "gencal(vis = '"+msName+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName+".tsys',\n"
    casaCmd = casaCmd + "  caltype = 'tsys')\n\n"

    calTableName1 = msName+'.tsys'

    chanEdge = 0.03125 # this is for 128ch/2GHz

    spwSpec = ''
    for i in sorted(tsysSpwInfo.keys()):
        if tsysSpwInfo[i]['numChans'] <= 256:
            if spwSpec != '': spwSpec = spwSpec+','
            chanEdge1 = chanEdge * tsysSpwInfo[i]['numChans'] / 128.
            spwSpec = spwSpec+str(i)+':0~'+str(np.longlong(tsysSpwInfo[i]['numChans']*chanEdge1-1))+';'+str(np.longlong(tsysSpwInfo[i]['numChans']-tsysSpwInfo[i]['numChans']*chanEdge1))+'~'+str(tsysSpwInfo[i]['numChans']-1)

    if spwSpec != '':
        casaCmd = casaCmd + "# Flagging edge channels\n\n"
        casaCmd = casaCmd + "flagdata(vis = '"+calTableName1+"',\n"
        casaCmd = casaCmd + "  mode = 'manual',\n"
        casaCmd = casaCmd + "  spw = '"+spwSpec+"',\n"
        casaCmd = casaCmd + "  flagbackup = False)\n\n"

    if isB2B:
        spwHighInt0 = sorted(sciSpwInfo.keys())[0]
        ref_freq = sciSpwInfo[spwHighInt0]['refFreq']
        if ref_freq > 600E9 and obsTimeStart < '2021-03-01T00:00:00': # Tsys mirroring was only needed up to March 2021
            casalog.post('B2B mode in Band 9 or 10 observed before March 2021: will add code to mirror Tsys', 'WARN')
            casaCmd = casaCmd + "# Mirroring Tsys\n"
            casaCmd = casaCmd + "tb.open('"+calTableName1+"')\n"
            casaCmd = casaCmd + "TsysSpws=np.unique(tb.getcol('SPECTRAL_WINDOW_ID'))\n"
            casaCmd = casaCmd + "tb.close()\n"
            casaCmd = casaCmd + "\n"  
            casaCmd = casaCmd + "TsysSpwSW=[]\n"
            casaCmd = casaCmd + "TsysSpwFreq=[]\n"
            casaCmd = casaCmd + "tb.open('"+calTableName1+"/SPECTRAL_WINDOW')\n"
            casaCmd = casaCmd + "for TsysSpw in TsysSpws:\n"
            casaCmd = casaCmd + "    spwFreqTsys = tb.getcell('REF_FREQUENCY',TsysSpw)\n"
            casaCmd = casaCmd + "    TsysSpwFreq.append(spwFreqTsys)\n"
            casaCmd = casaCmd + "    SWname = tb.getcell('NAME',TsysSpw)\n"
            casaCmd = casaCmd + "    TsysSpwSW.append(SWname[22:32])\n"
            casaCmd = casaCmd + "tb.close()\n"
            casaCmd = casaCmd + "\n"
            casaCmd = casaCmd + "ref_freq = "+str(ref_freq)+"\n"
            casaCmd = casaCmd + "\n"    
            casaCmd = casaCmd + "TsysIdsLF=np.where(np.array(TsysSpwFreq)<ref_freq/2)[0]\n"
            casaCmd = casaCmd + "TsysSpwsLF=TsysSpws[TsysIdsLF]\n"
            casaCmd = casaCmd + "TsysSpwsSWLF=np.array(TsysSpwSW)[TsysIdsLF]\n"
            casaCmd = casaCmd + "\n"
            casaCmd = casaCmd + "# now we want to create real and fake DSB pairings\n"
            casaCmd = casaCmd + "# e.g BB01 and SW01 with BB01 and SW02\n" 
            casaCmd = casaCmd + "baseBands = ['BB_1','BB_2','BB_3','BB_4']\n"
            casaCmd = casaCmd + "for baseBand in baseBands:\n"
            casaCmd = casaCmd + "    TsysPair=[spwNo for spwId,spwNo in enumerate(TsysSpwsLF) if baseBand in TsysSpwsSWLF[spwId]]\n"
            casaCmd = casaCmd + "    TsysData=[]\n"
            casaCmd = casaCmd + "    for TsysSpwPair in TsysPair:\n"
            casaCmd = casaCmd + "        tb.open('"+calTableName1+"', nomodify=False)\n"
            casaCmd = casaCmd + "        tsysdata = tb.query('SPECTRAL_WINDOW_ID == %s' %TsysSpwPair)\n"
            casaCmd = casaCmd + "        tsysdataSpw = tsysdata.getcol('FPARAM')\n"
            casaCmd = casaCmd + "        TsysData.append(tsysdataSpw)\n"
            casaCmd = casaCmd + "        tb.close()\n"
            casaCmd = casaCmd + "\n"
            casaCmd = casaCmd + "    if np.mean(TsysData[0][0][8:-8]) < np.mean(TsysData[1][0][8:-8]):\n"
            casaCmd = casaCmd + "        spwRep = 1\n"
            casaCmd = casaCmd + "        spwKeep = 0\n"
            casaCmd = casaCmd + "    else:\n"
            casaCmd = casaCmd + "        spwRep = 0\n"
            casaCmd = casaCmd + "        spwKeep =1\n"
            casaCmd = casaCmd + "    #now make the replacement\n"
            casaCmd = casaCmd + "    tb.open('"+calTableName1+"', nomodify=False)\n"
            casaCmd = casaCmd + "    tsysdataRep = tb.query('SPECTRAL_WINDOW_ID == %s' %TsysPair[spwRep])\n"
            casaCmd = casaCmd + "    tsysdataRep.putcol('FPARAM',TsysData[spwKeep])\n"
            casaCmd = casaCmd + "    tb.close()\n\n"


    if doplot:
        chanrange = '92.1875%'

        showimage = False
        for i in sorted(sciSpwInfo.keys()):
            if sciSpwInfo[i]['refFreq'] > 550e9: showimage = True

        if isB2B: # temporary workaround for a problem in plotbandpass
            casaCmd = casaCmd + "try: # protect against failure of plotbandpass in case of Tsys timestamp issues\n"
            casaCmd = casaCmd + "  if applyonly != True: aU.plotbandpass(caltable='%s', overlay='time', \n" %(calTableName1)
            casaCmd = casaCmd + "    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,\n"
            casaCmd = casaCmd + "    showatm=True,pwv='auto',chanrange='"+chanrange+"',showfdm=True, showBasebandNumber=True, showimage="+str(showimage)+", \n"
            casaCmd = casaCmd + "    field='', figfile='%s') \n\n" %(calTableName1+'.plots.overlayTime/'+calTableName1.split('/')[-1])
            casaCmd = casaCmd + "except:\n"
            casaCmd = casaCmd + "  print('Error in plotbandpass. Skipping the time-overlay Tsys plots.')\n"

        else:
            casaCmd = casaCmd + "if applyonly != True: aU.plotbandpass(caltable='%s', overlay='time', \n" %(calTableName1)
            casaCmd = casaCmd + "  xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,\n"
            casaCmd = casaCmd + "  showatm=True,pwv='auto',chanrange='"+chanrange+"',showfdm=True, showBasebandNumber=True, showimage="+str(showimage)+", \n"
            casaCmd = casaCmd + "  field='', figfile='%s') \n\n" %(calTableName1+'.plots.overlayTime/'+calTableName1.split('/')[-1])

    calTableName.append(calTableName1)
    if doplot:
        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+"', msName='"+msName+"', interactive=False) \n"

    return casaCmd

#################################

def doApplyAprioriCalTables(msName, tsys='', wvr='', antpos='', tsysmap='', tsysChanTol='', tsysPerField=False, valueMaps={}):
    """Generate code for the applycal step (apriori calibration: WVR, Tsys, and antpos) of a calibration script."""

    print('\n*** doApplyAprioriCalTables ***')
    print('Gathering information ...')

    casaCmd = ''

    if tsys=='' and wvr=='' and antpos=='': 
        casalog.post('ERROR: No cal table specified.', 'SEVERE')
        return False

    gainTable = []
    gainTable.append(tsys)
    gainTable.append(wvr)
    gainTable.append(antpos)
    gainTable = ['%s' %i for i in gainTable if i != '']

    spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS|CALIBRATE_DIFFGAIN', caching=True)
    spwIds = sorted(spwInfo.keys())
    spwIds1 = ','.join(['%s' %i for i in spwIds])

    if wvr=='':
        theinterp = "'linear,linear'" # default
    elif tsys=='' and antpos !='':
        theinterp = "['nearest, linear','']" # PRTSIR-16811: nearest interpol in time for WVR
    elif tsys!='' and antpos =='':
        theinterp = "['','nearest, linear']" # PRTSIR-16811: nearest interpol in time for WVR
    elif tsys=='' and antpos =='':
        theinterp = "['nearest, linear']" # PRTSIR-16811: nearest interpol in time for WVR
    else:
        theinterp = "['','nearest, linear','']" # PRTSIR-16811: nearest interpol in time for WVR

    if tsys=='':

        casaCmd = casaCmd + "\n\napplycal(vis = '"+msName+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwIds1+"',\n"
        casaCmd = casaCmd + "  gaintable = "+str(gainTable)+",\n"
        if re.search('^3.3', aU.getCasaVersion()) == None:
            casaCmd = casaCmd + "  interp = "+theinterp+",\n"
        else:
            casaCmd = casaCmd + "  interp = 'linear',\n"
        if aU.getCasaVersion() <= '4.2.2':
            casaCmd = casaCmd + "  calwt = False,\n"
        else:
            casaCmd = casaCmd + "  calwt = True,\n"
        casaCmd = casaCmd + "  flagbackup = False)\n"

    else:

        if re.search('^3.3', aU.getCasaVersion()) == None:
            if tsysmap == '':

                if tsysPerField == True:

                    casaCmd = casaCmd + "\n\nfrom almahelpers_localcopy import tsysspwmap\n"

                    if tsysChanTol == '':
                        casaCmd = casaCmd + "tsysmap = tsysspwmap(vis = '"+msName+"', tsystable = '"+tsys+"', perField=True)\n\n"
                    else:
                        tsysChanTol = int(tsysChanTol)
                        casaCmd = casaCmd + "tsysmap = tsysspwmap(vis = '"+msName+"', tsystable = '"+tsys+"', tsysChanTol = "+str(tsysChanTol)+", perField=True)\n\n"

                else:
                    if aU.getCasaVersion() < '5.9.9':
                        casaCmd = casaCmd + "\n\nfrom recipes.almahelpers import tsysspwmap\n"
                    else:
                        casaCmd = casaCmd + "\n\nfrom casarecipes.almahelpers import tsysspwmap\n"

                    if tsysChanTol == '':
                        casaCmd = casaCmd + "tsysmap = tsysspwmap(vis = '"+msName+"', tsystable = '"+tsys+"')\n\n"
                    else:
                        tsysChanTol = int(tsysChanTol)
                        casaCmd = casaCmd + "tsysmap = tsysspwmap(vis = '"+msName+"', tsystable = '"+tsys+"', tsysChanTol = "+str(tsysChanTol)+")\n\n"

            else:
                casaCmd = casaCmd + "\n\ntsysmap = "+str(tsysmap)+"\n\n"
        mytb = aU.createCasaTool(tbtool)
        mytb.open(msName+'/FIELD')
        sourceIds = mytb.getcol('SOURCE_ID')
        sourceNames = mytb.getcol('NAME')
        mytb.close()

        sourceIds1 = sorted(dict.fromkeys(sourceIds).keys())

        phaseCal = sfsdr.getPhaseCal(msName, valueMaps=valueMaps)

        intentSources = sfsdr.getIntentsAndSourceNames(msName)
        fieldIds = intentSources['CALIBRATE_ATMOSPHERE']['id']

        for i in sourceIds1:

            fieldIds1 = (np.where(sourceIds == i))[0]
            sourceName = sourceNames[fieldIds1[0]]

            if tsys != '':

                spwIds2 = []
                sourceIntents = []
                mymsmd = msmdtool()
                mymsmd.open(msName)
                for j in fieldIds1:
                    for k in mymsmd.spwsforfield(j):
                        spwIds2.append(k)
                    sourceIntents.append(mymsmd.intentsforfield(j))
                mymsmd.close()

                spwIds2 = np.unique(spwIds2)
                spwIds1 = ','.join(['%s' %j for j in spwIds if j in spwIds2])

                sourceIntents = np.unique(np.hstack(sourceIntents)).tolist()

                found = 0

                for j in range(len(sourceIntents)):
                    if re.search('CALIBRATE_ATMOSPHERE#[A-Z_]|CALIBRATE_WVR#[A-Z_]', sourceIntents[j]) == None:
                        found = 1
                        break

                if found == 0: continue

            if len(fieldIds1) > 1:
                j0 = 0
                fieldIds2 = str(fieldIds1[j0])
                for j in range(len(fieldIds1)-1):
                    if fieldIds1[j+1] == fieldIds1[j]+1: continue
                    fieldIds2 = fieldIds2 + '~' + str(fieldIds1[j])
                    j0 = j+1
                    fieldIds2 = fieldIds2 + ',' + str(fieldIds1[j0])
                fieldIds2 = fieldIds2 + '~' + str(fieldIds1[j+1])
            else:
                fieldIds2 = str(fieldIds1[0])

            fieldIds3 = [j for j in fieldIds1 if j in fieldIds]

            if len(fieldIds3) > 1:
                print('WARNING: Too many Tsys fields per source.')
                fieldIds3[0] = ','.join(['%d' %j for j in fieldIds3])

            if len(fieldIds3) == 0 or i in intentSources['OBSERVE_CHECK']['sourceid']:
                if i in intentSources['OBSERVE_TARGET']['sourceid']:

                    found = 0
                    for j in list(phaseCal.keys()):
                        if j == sourceName: continue
                        if phaseCal[j]['phaseCalId'] == phaseCal[sourceName]['phaseCalId']:
                            sourceIds2 = np.unique(phaseCal[j]['sciSourceIds'])
                            if len(sourceIds2) != 1: 
                                casalog.post('ERROR: no or more than one source ID','SEVERE')
                                return False
                            fieldIds3 = []
                            for k in range(len(intentSources['CALIBRATE_ATMOSPHERE']['sourceid'])):
                                if intentSources['CALIBRATE_ATMOSPHERE']['sourceid'][k] == sourceIds2:
                                    fieldIds3.append(intentSources['CALIBRATE_ATMOSPHERE']['id'][k])
                            if (len(fieldIds3) > 0):
                                fieldIds3[0] = ','.join(['%d' %k for k in fieldIds3])
                                found = 1
                                break
                    if found == 1:
                        casaCmd = casaCmd + "\n\n# Note: "+sourceName+" didn't have any Tsys measurement, so I used the one made on "+j+". This is probably Ok."
                    else:

                        if phaseCal[sourceName]['phaseCalId'] in fieldIds:
                            fieldIds3 = [phaseCal[sourceName]['phaseCalId']]
                            casaCmd = casaCmd + "\n\n# Note: "+sourceName+" didn't have any Tsys measurement, so I used the one made on "+phaseCal[sourceName]['phaseCalName']+". This is probably Ok."
                        else:
                            casaCmd = casaCmd + "\n\n# Warning: "+sourceName+" didn't have any Tsys measurement, and I couldn't find any close measurement. This is a science target, so this is probably *NOT* Ok."
                            continue

                elif i in intentSources['CALIBRATE_PHASE']['sourceid']:
                    if len(fieldIds1) != 1: 
                        casalog.post('ERROR: no or more than one field ID','SEVERE')
                        return False
                    found = 0
                    fieldIds3 = []
                    fieldIds3Names = []
                    for j in sorted(phaseCal.keys()):
                        if phaseCal[j]['phaseCalId'] == fieldIds1[0]:
                            sourceIds2 = np.unique(phaseCal[j]['sciSourceIds'])
                            if len(sourceIds2) != 1: 
                                casalog.post('ERROR: no or more than one source ID','SEVERE')
                                return False
                            for k in range(len(intentSources['CALIBRATE_ATMOSPHERE']['sourceid'])):
                                if intentSources['CALIBRATE_ATMOSPHERE']['sourceid'][k] == sourceIds2:
                                    fieldIds3.append(intentSources['CALIBRATE_ATMOSPHERE']['id'][k])
                                    fieldIds3Names.append(j)
                    if (len(fieldIds3) > 0):
                        fieldIds3[0] = ','.join(['%d' %k for k in sorted(fieldIds3)])
                        fieldIds3Names = ','.join(fieldIds3Names)
                        found = 1


                    if found == 1:
                        casaCmd = casaCmd + "\n\n# Note: "+sourceName+" didn't have any Tsys measurement, so I used the one made on "+fieldIds3Names+". This is probably Ok."
                    else:
                        casaCmd = casaCmd + "\n\n# Warning: "+sourceName+" didn't have any Tsys measurement, and I couldn't find any close measurement. This is a phase calibrator, so this is probably *NOT* Ok."
                        continue
                elif i in intentSources['OBSERVE_CHECK']['sourceid']:
                    if len(fieldIds1) != 1: 
                        casalog.post('ERROR: no or more than one field ID','SEVERE')
                        return False
                    found = 0
                    for j in list(phaseCal.keys()):
                        if j == sourceName: continue
                        if phaseCal[j]['phaseCalId'] == phaseCal[sourceName]['phaseCalId']:
                            found = 1
                            break
                    if found == 1:
                        sourceIds2 = np.unique(phaseCal[j]['sciSourceIds'])
                        if len(sourceIds2) != 1: 
                            casalog.post('ERROR: no or more than one source ID','SEVERE')
                            return False
                        fieldIds3 = []
                        for k in range(len(intentSources['CALIBRATE_ATMOSPHERE']['sourceid'])):
                            if intentSources['CALIBRATE_ATMOSPHERE']['sourceid'][k] == sourceIds2:
                                fieldIds3.append(intentSources['CALIBRATE_ATMOSPHERE']['id'][k])

                        if (len(fieldIds3) > 0):
                            fieldIds3[0] = ','.join(['%d' %k for k in fieldIds3])
                        casaCmd = casaCmd + "\n\n# Note: "+sourceName+" didn't have any Tsys measurement, so I used the one made on "+j+". This is probably Ok."
                    else:
                        casaCmd = casaCmd + "\n\n# Warning: "+sourceName+" didn't have any Tsys measurement, and I couldn't find any close measurement. This is a check source, so this is probably *NOT* Ok."
                        continue
                else:

                    found = 0

                    if len(fieldIds1) == 1:
                        fieldIds1 = fieldIds1[0]
                        if phaseCal is not None:
                            for j in phaseCal:
                                if phaseCal[j]['phaseCalId'] == fieldIds1:
                                    sciSourceIds1 = phaseCal[j]['sciSourceIds']
                                    sciSourceIds1 = sorted(dict.fromkeys(sciSourceIds1).keys())
                                    if len(sciSourceIds1) == 1:
                                        sciSourceIds1 = sciSourceIds1[0]
                                        ij = np.where(np.array(intentSources['CALIBRATE_ATMOSPHERE']['sourceid']) == sciSourceIds1)[0]
                                        if len(ij) == 1:
                                            ij = ij[0]
                                            fieldIds3 = [intentSources['CALIBRATE_ATMOSPHERE']['id'][ij]]
                                            found = 1
                                            break

                    if found == 0:

                        if tsys != '':

                            fieldIds4 = (np.where(sourceIds == i))[0]
                            mymsmd = msmdtool()
                            mymsmd.open(msName)
                            fieldIds5 = mymsmd.fieldsforname(sourceName)

                            fieldIds5 = [j for j in fieldIds5 if j not in fieldIds4]

                            fieldIds6 = []
                            for j in fieldIds5:
                                fieldIntents = mymsmd.intentsforfield(j)
                                if 'CALIBRATE_ATMOSPHERE#ON_SOURCE' in fieldIntents and 'CALIBRATE_ATMOSPHERE#OFF_SOURCE' in fieldIntents:
                                    fieldIds6.append(j)

                            mymsmd.close()

                            if len(fieldIds6) != 0:

                                fieldIds3 = [','.join(['%d' %k for k in fieldIds6])]
                                casaCmd = casaCmd + "\n\n# Note: "+sourceName+" didn't have any Tsys measurement, so I used the one made on another source with the same name. This is probably Ok."
                                found = 1

                        if found == 0:

                            casaCmd = casaCmd + "\n\n# Note: "+sourceName+" didn't have any Tsys measurement, and I couldn't find any close measurement. But this is not a science target, so this is probably Ok."
                            continue

            gainField = []
            for j in range(len(gainTable)): gainField.append('')

            if (len(fieldIds3) > 0):
                gainField[0] = str(fieldIds3[0])
            else:
                gainField[0] = intentSources['CALIBRATE_ATMOSPHERE']['idstring'][0]

            gainSpwMap = []
            for j in range(len(gainTable)): gainSpwMap.append("[]")

            if tsysPerField == True:
                gainSpwMap[0] = "tsysmap['"+str(gainField[0])+"']"
            else:
                gainSpwMap[0] = 'tsysmap'

            gainSpwMap = ','.join(gainSpwMap)

            casaCmd = casaCmd + "\n\napplycal(vis = '"+msName+"',\n"
            casaCmd = casaCmd + "  field = '"+fieldIds2+"',\n"
            casaCmd = casaCmd + "  spw = '"+spwIds1+"',\n"
            casaCmd = casaCmd + "  gaintable = "+str(gainTable)+",\n"
            casaCmd = casaCmd + "  gainfield = "+str(gainField)+",\n"
            if re.search('^3.3', aU.getCasaVersion()) == None:
                casaCmd = casaCmd + "  interp = "+theinterp+",\n"
            else:
                casaCmd = casaCmd + "  interp = 'linear',\n"
            if re.search('^3.3', aU.getCasaVersion()) == None: casaCmd = casaCmd + "  spwmap = ["+gainSpwMap+"],\n"
            casaCmd = casaCmd + "  calwt = True,\n"
            casaCmd = casaCmd + "  flagbackup = False)\n\n"

        casaCmd = casaCmd + "\n\nif applyonly != True: es.getCalWeightStats('"+msName+"', 0.2) \n"

    return casaCmd


##################################

def doSplitOut(msName, msName1='', outMsName='', splitMyScienceSpw=False, timebin=0., iHaveSplitMyScienceSpw=False, allowHybrid=True, intentsToDiscard='', reindexMyScienceSpw=False):
    """Generate code for the split-out the corrected data."""

    print('\n*** doSplitOut ***')
    print('Gathering information ...')

    if aU.getCasaVersion() >= '5.4':
        useMStransform = True
    else:
        useMStransform = False

    casaCmd = ''

    if msName1 == '': msName1 = msName
    if outMsName == '': outMsName = msName1+'.split'

    spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS|CALIBRATE_DIFFGAIN', caching=True)
    spwIds = sorted(spwInfo.keys())

    if intentsToDiscard != '':

        mymsmd = msmdtool()
        mymsmd.open(msName)

        if iHaveSplitMyScienceSpw == True:

            listOfIntents = []

            for i in spwIds:
                for j in mymsmd.intentsforspw(i):
                    listOfIntents.append(j)

            listOfIntents = np.unique(listOfIntents).tolist()

        else:

            listOfIntents = mymsmd.intents()

        mymsmd.close()

        listOfIntents = [i for i in listOfIntents if re.search(intentsToDiscard, i, re.IGNORECASE) == None]

    if iHaveSplitMyScienceSpw == True: 
        spwIds = list(range(len(spwIds)))
    spwIds = ['%d' %i for i in spwIds]
    spwIds = ','.join(spwIds)

    numChans = []
    for i in sorted(spwInfo.keys()): numChans.append(spwInfo[i]['numChans'])
    if max(numChans) <= 256 and timebin != 0: casaCmd = casaCmd + "# Important note: the correlator mode for this dataset was TDM, you may want not to do any time averaging.\n\n"

    if allowHybrid != True:
        mytb = aU.createCasaTool(tbtool)
        mytb.open(msName+'/ANTENNA')
        antDiam = mytb.getcol('DISH_DIAMETER')
        antNames = mytb.getcol('NAME')
        mytb.close()

        antDiam1 = np.unique(antDiam)

        antNames2 = ''

        if len(antDiam1) != 1:

            casalog.post('This dataset has mixed antenna diameters (hybrid) and you set allowHybrid=False .','WARN')
            casalog.post('Will determine the dominant antenna diameter and only split out those antennas.','WARN')

            numAnts = {}
            for i in range(len(antDiam1)):
                numAnts[antDiam1[i]] = len(np.where(antDiam == antDiam1[i])[0])

            numAntsMax = max(numAnts.values())

            antDiam2 = []
            for i in range(len(antDiam1)):
                if numAnts[antDiam1[i]] == numAntsMax:
                    antDiam2.append(antDiam1[i])
            if len(antDiam2) != 1: 
                casalog.post('ERROR: There are equally many antennas with diameters '+str(antDiam2)+'. Cannot decide on dominant type.','SEVERE')
                return False
            antDiam2 = antDiam2[0]

            ij = np.where(antDiam == antDiam2)
            antNames1 = antNames[ij]

            antNames2 = []
            for i in range(len(antNames1)):
                antNames3 = re.findall('^[a-z]+', antNames1[i], re.IGNORECASE)
                if len(antNames3) == 0: continue
                antNames2.append(antNames3[0]+'*')
            antNames2 = ','.join(np.unique(antNames2).tolist())+'&'

    casaCmd = casaCmd + "os.system('rm -rf %s') \n"%(outMsName)
    casaCmd = casaCmd + "os.system('rm -rf %s.flagversions') \n\n"%(outMsName)
    if intentsToDiscard != '': casaCmd = casaCmd + "listOfIntents = %s\n\n"%(pprint.pformat(listOfIntents))
    if useMStransform == True:
        casaCmd = casaCmd + "mstransform(vis = '"+msName1+"',\n"
    else:
        casaCmd = casaCmd + "split(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  outputvis = '"+outMsName+"',\n"
    casaCmd = casaCmd + "  datacolumn = 'corrected',\n"
    if splitMyScienceSpw == True:
        casaCmd = casaCmd + "  spw = '"+spwIds+"',\n"
        if useMStransform == True: casaCmd = casaCmd + "  reindex = "+str(reindexMyScienceSpw)+",\n"
    if intentsToDiscard != '': casaCmd = casaCmd + "  intent = ','.join(listOfIntents),\n"        
    if timebin != 0.:
        if useMStransform == True: casaCmd = casaCmd + "  timeaverage = True,\n"
        casaCmd = casaCmd + "  timebin = '"+str(timebin)+"s',\n"
    if allowHybrid != True and antNames2 != '': 
        casaCmd = casaCmd + "  antenna = '"+antNames2+"',\n"
        casalog.post('Split-out step will only write antennas with diameter '+str(antDiam2)+": "+antNames2,'WARN')

    casaCmd = casaCmd + "  keepflags = True)\n\n"

    return casaCmd


###################################

def doSaveFlags(msName, name=''):
    """Generate code for the flag-saving step of a calibration script in preparation of
    a subsequent step which modifies the flags.

    name - The flag version name (obligatory)"""

    print('\n*** doSaveFlags ***')
    print('Gathering information ...')

    if name == '': 
        casalog.post('ERROR: Missing version name.','SEVERE')
        return False

    casaCmd = ''

    if name == 'Original':
        casaCmd = casaCmd + "\nif not os.path.exists('"+msName+".flagversions/Original.flags'):\n"
        casaCmd = casaCmd + "  flagmanager(vis = '"+msName+"',\n"
        casaCmd = casaCmd + "    mode = 'save',\n"
        casaCmd = casaCmd + "    versionname = '"+name+"')\n\n"
    else:
        casaCmd = casaCmd + "\nflagmanager(vis = '"+msName+"',\n"
        casaCmd = casaCmd + "  mode = 'save',\n"
        casaCmd = casaCmd + "  versionname = '"+name+"')\n\n"

    return casaCmd

####################################

def doClearPointingTable(msName):

    print('\n*** doClearPointingTable ***')

    casaCmd = ''

    casaCmd = casaCmd + "tb.open('"+msName+"/POINTING', nomodify = False)\n"
    casaCmd = casaCmd + "a = tb.rownumbers()\n"
    casaCmd = casaCmd + "tb.removerows(a)\n"
    casaCmd = casaCmd + "tb.close()\n"

    return casaCmd

####################################

def doInitialFlagging(msName, msName1='', chanEdge=0.0625, thresh=0.2, iHaveSplitMyScienceSpw=False):
    """Generate code for the initial flagging step of a calibration script."""

    print('\n*** doInitialFlagging ***')
    print('Gathering information ...')

    specLines = {'Neptune': [[114.00,116.50], [227.00,234.50], [340.00,351.50], [455.00,467.50], [686.00,696.50], [803.00,810.50]], # CO
        'Titan': [[114.93,115.66], [229.51,231.71], [343.86,347.58], [458.34,463.74], [687.83,694.58], [803.55,809.76], # CO
        [110.19,110.21], [220.30,220.50], [330.39,330.78], [440.47,441.10], [660.68,661.46], [770.74,771.63], [880.83,881.72], # 13CO
        [88.46,88.80], [176.75,177.78], [264.99,266.78], [353.33,355.68], [441.79,444.47], [618.52,622.09], [707.09,710.66], [795.65,799.22], [883.92,887.75], # HCN
        [86.05,86.06], [172.07,172.14], [258.07,258.24], [430.07,430.40], [602.02,602.53], [773.97,774.55], [859.93,860.52], # HC15N
        [86.34,86.34], [172.65,172.71], [258.94,259.09], [431.49,431.83], [604.04,604.49], [776.56,777.08], [862.81,863.33] # H13CN
        ]}

    if msName1 == '': msName1 = msName

    casaCmd = ''

    casaCmd = casaCmd + "# Flagging shadowed data\n\n"
    casaCmd = casaCmd + "flagdata(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  mode = 'shadow',\n"
    casaCmd = casaCmd + "  flagbackup = False)\n\n"

    spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS|CALIBRATE_DIFFGAIN', caching=True)
    spwIds = sorted(spwInfo.keys())

    sciNumChans = []
    for i in sorted(spwInfo.keys()): sciNumChans.append(spwInfo[i]['numChans'])
    if len(list(dict.fromkeys(sciNumChans).keys())) != 1:
        print('WARNING: This seems to be a mixed-mode dataset. Please check the script carefully.')

    spwSpec = ''

    mymsmd = msmdtool()
    mymsmd.open(msName)
    TDMspws = mymsmd.tdmspws()
    mymsmd.close()

    for i in range(len(spwIds)):
        if spwIds[i] in TDMspws:
            if spwSpec != '': spwSpec = spwSpec+','
            if iHaveSplitMyScienceSpw == True:
                spwSpec = spwSpec+str(i)
            else:
                spwSpec = spwSpec+str(spwIds[i])
            spwSpec = spwSpec+':0~'+str(np.longlong(sciNumChans[i]*chanEdge-1))+';'+str(np.longlong(sciNumChans[i]-sciNumChans[i]*chanEdge))+'~'+str(sciNumChans[i]-1)


    if spwSpec != '':
        casaCmd = casaCmd + "# Flagging edge channels\n\n"
        casaCmd = casaCmd + "flagdata(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  mode = 'manual',\n"
        casaCmd = casaCmd + "  spw = '"+spwSpec+"',\n"
        casaCmd = casaCmd + "  flagbackup = False)\n\n"

    spwInfo1 = sfsdr.getSpwInfo(msName, intent='CALIBRATE_AMPLI|CALIBRATE_FLUX', caching=True)

    if len(spwInfo1) != 0:

        spwIds = sorted(spwInfo1.keys())

        intentSources = sfsdr.getIntentsAndSourceNames(msName)
        calFieldIds = intentSources['CALIBRATE_AMPLI']['id'] + intentSources['CALIBRATE_FLUX']['id']
        calFieldIds = [i for i in calFieldIds if i != '']
        calFieldNames = intentSources['CALIBRATE_AMPLI']['name'] + intentSources['CALIBRATE_FLUX']['name']
        calFieldNames = [i for i in calFieldNames if i != '']

        for i in range(len(calFieldIds)):
            if calFieldNames[i] in list(specLines.keys()):
                for j in specLines[calFieldNames[i]]:
                    for k in range(len(spwIds)):
                        spwId1 = int(sorted(spwInfo1.keys())[k])
                        chanRange = aU.getChanRangeFromFreqRange(vis = msName, spwid = spwId1, minf = j[0]*1.e9, maxf = j[1]*1.e9)
                        if chanRange == [-1, -1]: continue

                        if iHaveSplitMyScienceSpw == True:
                            spwId2 = sorted(spwInfo.keys()).index(spwId1)
                        else:
                            spwId2 = spwId1

                        if (chanRange[1]-chanRange[0]) / (spwInfo[spwId1]['numChans']*1.) > thresh: print('# Warning: more than '+str(thresh*100)+'% of spw '+str(spwId2)+' on '+calFieldNames[i]+' will be flagged due to atmospheric line.')
                        spwSpec = str(spwId2)+':'+str(chanRange[0])+'~'+str(chanRange[1])

                        casaCmd = casaCmd + "# Flagging atmospheric line(s)\n\n"
                        casaCmd = casaCmd + "flagdata(vis = '"+msName1+"',\n"
                        casaCmd = casaCmd + "  mode = 'manual',\n"
                        casaCmd = casaCmd + "  spw = '"+spwSpec+"',\n"
                        casaCmd = casaCmd + "  field = '"+str(calFieldIds[i])+"',\n"
                        casaCmd = casaCmd + "  flagbackup = False)\n\n"

    return casaCmd

#####################################

def doRunSetjy(msName, msName1='', iHaveSplitMyScienceSpw=False, useCalibratorService=False, calibratorServiceURL=None, 
               isB2B=False, valueMaps={}):
    """Generate code for the setjy step of a calibration script,
    i.e. for setting a model for the flux calibrator(s).

    Return the generated code and the list of the name(s) of the selected flux calibrators

    """

    print('\n*** doRunSetjy ***')
    print('Gathering information ...')

    if msName1 == '': msName1 = msName

    casaCmd = ''
    theFluxCalNames = [] # store the names of the fields which are actually used as fluxcal (in order of field ID)

    print('Decide on model or quasar fluxcal ...')

    fieldIds = sfsdr.getFieldsForSetjy(msName)
    mytb = aU.createCasaTool(tbtool)

    mytb.open(msName+'/FIELD')
    fieldNames = mytb.getcol('NAME')
    ephemerisIds = mytb.getcol('EPHEMERIS_ID')
    phaseDirKeywords = mytb.getcolkeywords('PHASE_DIR')
    phaseDirRef = mytb.getcol('PhaseDir_Ref')
    mytb.close()

    ij = np.where(phaseDirKeywords['MEASINFO']['TabRefTypes'] == 'ICRS')[0][0]
    icrscode = phaseDirKeywords['MEASINFO']['TabRefCodes'][ij]

    if fieldIds != []: # there are non-quasar flux calibrators

        print('Model fluxcal ...')

        fieldNames1 = ['%s' %fieldNames[i] for i in fieldIds] # the names of the non-quasar flux calibrators
        fieldNames = ','.join(fieldNames1)
        fieldIds = ['%s' %i for i in fieldIds]
        fieldIds1 = ','.join(fieldIds)

        spwInfo = sfsdr.getSpwInfo(msName, intent='CALIBRATE_AMPLI|CALIBRATE_FLUX', caching=True)
        spwIds = sorted(spwInfo.keys())

        if iHaveSplitMyScienceSpw == True:
            spwInfo1 = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS|CALIBRATE_DIFFGAIN', caching=True)
            spwIds1 = sorted(spwInfo1.keys())
            spwIds = [spwIds1.index(i) for i in spwIds]
            
        spwIds = ['%s' %i for i in spwIds]
        spwIds = ','.join(spwIds)

        if aU.getCasaVersion() == '4.5.0':

            for i in range(len(fieldIds)):

                casaCmd = casaCmd + "setjy(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  field = '"+fieldIds[i]+"', # "+fieldNames1[i]+"\n"
                theFluxCalNames.append(fieldNames1[i])
                casaCmd = casaCmd + "  spw = '"+spwIds+"',\n"

                if ephemerisIds[int(fieldIds[i])] != -1 and phaseDirRef[int(fieldIds[i])] == icrscode:
                    casaCmd = casaCmd + "  standard = 'Butler-JPL-Horizons 2012',\n"
                    casaCmd = casaCmd + "  useephemdir = True)\n\n"
                else:
                    casaCmd = casaCmd + "  standard = 'Butler-JPL-Horizons 2012')\n\n"

        else:

            casaCmd = casaCmd + "setjy(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  field = '"+fieldIds1+"', # "+fieldNames+"\n"
            theFluxCalNames = fieldNames1.copy()
            casaCmd = casaCmd + "  spw = '"+spwIds+"',\n"
            haveEphemFluxcal=False
            for i in range(len(fieldIds)):
                if ephemerisIds[int(fieldIds[i])] != -1 and phaseDirRef[int(fieldIds[i])] == icrscode:
                    haveEphemFluxcal=True
                    casalog.post('There are ephemeris objects among the flux calibrators: field id '+str(fieldIds[i]), 'WARN')

            casaCmd = casaCmd + "  usescratch = True,\n" # because of bug in virt. model (see SCIREQ-2035)
            casaCmd = casaCmd + "  standard = 'Butler-JPL-Horizons 2012')\n\n"

        casaCmd = casaCmd + "if applyonly != True:\n"
        casaCmd = casaCmd + "  os.system('rm -rf %s.setjy.field*.png') \n"%(msName1)
        casaCmd = casaCmd + "  for i in "+str(fieldIds)+":\n"
        casaCmd = casaCmd + "    plotms(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "      xaxis = 'uvdist',\n"
        casaCmd = casaCmd + "      yaxis = 'amp',\n"
        casaCmd = casaCmd + "      ydatacolumn = 'model',\n"
        casaCmd = casaCmd + "      field = str(i),\n"
        casaCmd = casaCmd + "      spw = '"+spwIds+"',\n"
        casaCmd = casaCmd + "      avgchannel = '9999',\n"
        casaCmd = casaCmd + "      coloraxis = 'spw',\n"
        casaCmd = casaCmd + "      plotfile = '"+msName1+".setjy.field'+i+'.png')\n"

    else: # fieldIds == [], i.e. we have a quasar fluxcal

        print('Quasar fluxcal ...')

        intentSources = sfsdr.getIntentsAndSourceNames(msName)
        myIntentSources = [intentSources['CALIBRATE_FLUX']]

        explicitSpw = [[],[]]
        if isB2B:
            myIntentSources.append(intentSources['CALIBRATE_DIFFGAIN'])
            allDGSpws = myIntentSources[1]['spw']
            spwInfo1 = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS|CALIBRATE_DIFFGAIN', caching=True)
            spwIds1 = sorted(spwInfo1.keys())

            print("fluxcal SPWs:  ", myIntentSources[0]['spw'])
            print("diffgain SPWs: ", allDGSpws)

            if myIntentSources[0]['spw'] == allDGSpws: # fluxcal covers HF and LF
                myIntentSources = [intentSources['CALIBRATE_FLUX']]
                explicitSpw = [[]]
                for i in allDGSpws:
                    if iHaveSplitMyScienceSpw == True:
                        if i in spwIds1:
                            explicitSpw[0].append(spwIds1.index(i))
                    else:
                        if i in spwIds1:
                            explicitSpw[0].append(i)
                print('B2B: Will use field id '+str(myIntentSources[0]['sourceid'])+' as flux cal for all science and diffgain SPWs: '+str(explicitSpw[0]))
            else:
                for i in allDGSpws:
                    if iHaveSplitMyScienceSpw == True:
                        if  i in myIntentSources[0]['spw'] and i in spwIds1:
                            explicitSpw[0].append(spwIds1.index(i))
                        elif i in spwIds1:
                            explicitSpw[1].append(spwIds1.index(i))
                    else:
                        if  i in myIntentSources[0]['spw'] and i in spwIds1:
                            explicitSpw[0].append(i)
                        elif i in spwIds1:
                            explicitSpw[1].append(i)
                        
                for j in range(0,2):
                    print('B2B: Will use field id '+str(myIntentSources[j]['sourceid'])+' as flux cal for SPWs '+str(explicitSpw[j]))

        else: # not B2B
            print("fluxcal SPWs:  ", myIntentSources[0]['spw'])
            spwInfo1 = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
            spwIds1 = sorted(spwInfo1.keys())

            

        for ii in range(len(myIntentSources)):

            if len(myIntentSources[ii]['sourceid']) != 0:

                fluxCalSourceId = myIntentSources[ii]['sourceid']
                fluxCalSourceNames = myIntentSources[ii]['name']

                found = 0
                
                sourceFluxesPerSPW = {}
                confimedFluxCalNames = set([])
                universal_spix = True
                
                if useCalibratorService:
                    print('Running aU.getALMAFluxForMS per science SPW using the calibrator service ...')

                    for myspw in myIntentSources[0]['spw']:
                        if myspw in spwIds1:
                            print('   SPW ', myspw)
                            sFluxes = aU.getALMAFluxForMS(msName, useCalibratorService=useCalibratorService, calibratorServiceURL=calibratorServiceURL, spw=myspw)
                            if len(sFluxes) != 0:
                                for j in list(sFluxes.keys()):
                                    if j not in fluxCalSourceNames:
                                        sFluxes.pop(j)
                                fluxCalIds = []
                                for j in sFluxes.keys():
                                    fluxCalIds.append(list(fieldNames).index(j))
                                fluxCalId = sorted(fluxCalIds)[0] # pick the fluxcal with the lowest ID
                                fluxCalSourceName = fieldNames[fluxCalId]
                                confimedFluxCalNames.add(fluxCalSourceName)
                                
                                if len(sFluxes.keys()) > 1:
                                    casalog.post("THERE IS MORE THAN ONE FLUX CALIBRATOR. WILL PICK THE FIRST ONE: "+fluxCalSourceName+". THIS MAY BE WRONG.", 'WARN')
                                for j in list(sFluxes.keys()):
                                    if j != fluxCalSourceName:
                                        sFluxes.pop(j)
                            else: # sFluxes empty
                                casalog.post("   ERROR: There is no usable flux catalog information. Will try to continue ...",'WARN')
                                        
                            sourceFluxesPerSPW[myspw] = sFluxes

                    if len(confimedFluxCalNames) == 0:
                        casalog.post("ERROR: There is no usable flux calibrator.",'SEVERE')
                        return False
                    if len(confimedFluxCalNames) > 1:
                        casalog.post("ERROR: (internal) a different flux calibrator was choses for different SPWs.",'SEVERE')
                        return False
                    
                            
                    # check if spix values are the same for all SPWs
                    myFluxCalSourceName = list(confimedFluxCalNames)[0]
                    print('\n*** Flux calibrator '+myFluxCalSourceName+':')
                    firstspw = list(sourceFluxesPerSPW.keys())[0]
                    firstspix = round(sourceFluxesPerSPW[firstspw][myFluxCalSourceName]['spectralIndex'],6)
                    for myspw in sourceFluxesPerSPW.keys():
                        myspix = round(sourceFluxesPerSPW[myspw][myFluxCalSourceName]['spectralIndex'],6)
                        print('     SPW: '+str(myspw)+', Spix: '+str(myspix))
                        if not myspix == firstspix:
                            universal_spix = False
                                   
                    if universal_spix:
                        print('*** All SPWs share the same spix. Re-querying calibrator service for entire spectral range ...')
                        sFluxes = aU.getALMAFluxForMS(msName, useCalibratorService=useCalibratorService, calibratorServiceURL=calibratorServiceURL)
                        if len(sFluxes) != 0:
                            sourceFluxesPerSPW = {}
                            sourceFluxesPerSPW['allspws'] = sFluxes
                            
                else:
                    print('Running aU.getALMAFluxForMS without using the calibrator service ...')
                        
                    sourceFluxesPerSPW['allspws'] = aU.getALMAFluxForMS(msName)
                    
                    
                for myspw in sourceFluxesPerSPW.keys(): # loop over the getALMAFluxForMS return values for each science SPW
                    if myspw != 'allspws' and myspw not in spwIds1:
                        continue   # only do the science SPWs
                        
                    sourceFluxes = sourceFluxesPerSPW[myspw]
                    
                    if len(sourceFluxes) != 0:

                        for j in list(sourceFluxes.keys()):
                            if j not in fluxCalSourceNames:
                                sourceFluxes.pop(j)

                        if len(sourceFluxes) != 0:

                            found = 1
                            fluxCalSourceNames = list(sourceFluxes.keys())
                            fluxCalSourceName = fluxCalSourceNames[0]
                            
                            if len(fluxCalSourceNames) > 1:
                                casalog.post("THERE IS MORE THAN ONE FLUX CALIBRATOR. WILL PICK THE FIRST ONE: "+fluxCalSourceName+". THIS MAY BE WRONG.", 'WARN')
                                
                            fluxCalId = list(fieldNames).index(fluxCalSourceName)
                            haveEphemFluxcal=False
                            if ephemerisIds[fluxCalId] != -1 and phaseDirRef[fluxCalId] == icrscode:
                                haveEphemFluxcal=True
                                casalog.post('  The flux calibrator is an ephemeris object: '+fluxCalSourceName, 'WARN')
                            elif 'DataConditions' in sourceFluxes[fluxCalSourceName].keys(): # aU supports data conditions check
                                for mySName in fluxCalSourceNames:
                                    myCond = str(sourceFluxes[mySName]['DataConditions'])
                                    casalog.post("   Flux Calibrator "+mySName+" has condition "+myCond, 'INFO')

                                fluxCalConditions = str(sourceFluxes[fluxCalSourceName]['DataConditions'])
                                # value[0]: number of available measurements used, where 9 means 9 or more
                                # value[1]: 1 if measurements exist from at least two distinct bands, 0 otherwise
                                # value[2]: 1 if there is at least one measurement on either side of the selected date, 0 otherwise)
                                if len(fluxCalConditions)>=3:
                                    casaCmd = casaCmd + "\n# URL for catalog access: "+sourceFluxes[fluxCalSourceName]['url']+"\n"
                                    casaCmd = casaCmd + "# Number of available measurements used (where 9 means 9 or more): "+str(fluxCalConditions[0])+"\n"
                                    casaCmd = casaCmd + "# Measurements were available in more than one band: "+str(fluxCalConditions[1]=='1')+"\n"
                                    casaCmd = casaCmd + "# Measurements bracketed the observation date: "+str(fluxCalConditions[2]=='1')+"\n\n"
                                    if myspw == 'allspws': # setjy for all SPWs at once
                                        casaCmd = casaCmd +"# The catalog value for the flux calibrator spectral index was constant over all SPWs.\n\n"
                                else:
                                    casalog.post('Invalid condition code returned by aU.getALMAFluxForMS for flux calibrator '+fluxCalSourceName+': '+fluxCalConditions, 'WARN')

                            casaCmd = casaCmd + "setjy(vis = '"+msName1+"',\n"
                            casaCmd = casaCmd + "  standard = 'manual',\n"
                            casaCmd = casaCmd + "  field = '"+fluxCalSourceName+"',\n"
                            if fluxCalSourceName not in theFluxCalNames:
                                theFluxCalNames.append(fluxCalSourceName)

                            if myspw == 'allspws': # setjy for all SPWs at once
                                if isB2B and explicitSpw[ii] != []:
                                    casaCmd = casaCmd + "  spw = '"+','.join(str(n) for n in explicitSpw[ii])+"',\n"
                            else:
                                if iHaveSplitMyScienceSpw == True:                                    
                                    casaCmd = casaCmd + "  spw = '"+str(spwIds1.index(myspw))+"',\n"
                                else:
                                    casaCmd = casaCmd + "  spw = '"+str(myspw)+"',\n"

                            casaCmd = casaCmd + "  usescratch = True,\n" # because of bug in virt. model (see SCIREQ-2035)

                            casaCmd = casaCmd + "  fluxdensity = ["+str(sourceFluxes[fluxCalSourceName]['fluxDensity'])+", 0, 0, 0],\n"
                            casaCmd = casaCmd + "  spix = "+str(round(sourceFluxes[fluxCalSourceName]['spectralIndex'],6))+",\n"
                            casaCmd = casaCmd + "  reffreq = '"+str(sourceFluxes[fluxCalSourceName]['frequency']/1.e9)+"GHz')\n\n"

                            for tag in ['fluxDensityUncertainty', 'meanAge']: casaCmd = casaCmd + "# "+tag+" = "+str(sourceFluxes[fluxCalSourceName][tag])+"\n"

                    if found == 0:

                        fluxCalSourceId = myIntentSources[ii]['sourceid']

                        sourceFluxes = sfsdr.getFluxesFromSourceTable(msName)

                        if len(fluxCalSourceId) > 1:
                            print("WARNING: THERE IS MORE THAN ONE FLUX CALIBRATOR. I WILL PICK THE FIRST ONE. THIS MAY BE WRONG.")
                            fluxCalSourceId = [j for j in fluxCalSourceId if j in list(sourceFluxes.keys())]

                        if len(fluxCalSourceId) == 0: 
                            casalog.post("ERROR: There is no flux calibrator.",'SEVERE')
                            return False

                        fluxCalSourceId = fluxCalSourceId[0]

                        fluxCalSourceName = myIntentSources[ii]['name'][myIntentSources[ii]['sourceid'].index(fluxCalSourceId)]

                        if len(myIntentSources[ii]['sourceid']) > 1:
                            mytb.open(msName+'/FIELD')
                            tb1 = mytb.query('SOURCE_ID == '+str(fluxCalSourceId))
                            fluxCalFieldIds1 = tb1.rownumbers().tolist()
                            tb1.close()
                            mytb.close()
                            fluxCalFieldIds = [j for j in fluxCalFieldIds1 if j in myIntentSources[i]['id']]
                        else:
                            fluxCalFieldIds = myIntentSources[ii]['id']

                        if fluxCalSourceId in sourceFluxes:

                            if fluxCalSourceName != sourceFluxes[fluxCalSourceId]['sourceName']: 
                                casalog.post("ERROR: Source names do not match.",'SEVERE')
                                return False

                            haveEphemFluxcal=False
                            for k in range(len(fluxCalFieldIds)):
                                if ephemerisIds[int(fluxCalFieldIds[k])] != -1 and phaseDirRef[int(fluxCalFieldIds[k])] == icrscode:
                                    haveEphemFluxcal=True
                                    casalog.post('There are ephemeris objects among the flux calibrators: field id '+str(fluxCalFieldIds[k]), 'WARN')

                            fluxCalFieldIds = ['%s' %k for k in fluxCalFieldIds]
                            fluxCalFieldIds = ','.join(fluxCalFieldIds)

                            if msName in valueMaps.keys():
                                vm = valueMaps[msName]
                                print('Using canned ValueMap.')
                            else:
                                vm = aU.ValueMapping(msName)
                                valueMaps[msName] = vm

                            spwInfo = sfsdr.getSpwInfo(msName, intent='CALIBRATE_FLUX', caching=True)
                            spwIds = sorted(spwInfo.keys())

                            spwMeanFreq = []
                            for j in spwIds: spwMeanFreq.append(vm.spwInfo[j]['meanFreq'])

                            if iHaveSplitMyScienceSpw == True: 
                                spwIds = list(range(len(spwIds)))

                            for j in range(len(spwIds)):

                                frequency1 = []
                                for k in sourceFluxes[fluxCalSourceId]['frequency']: frequency1.append(abs(k-spwMeanFreq[j]))
                                ij = frequency1.index(min(frequency1))

                                frequency1 = sourceFluxes[fluxCalSourceId]['frequency'][ij]
                                flux1 = sourceFluxes[fluxCalSourceId]['flux'][ij]

                                casaCmd = casaCmd + "setjy(vis = '"+msName1+"',\n"
                                casaCmd = casaCmd + "  field = '"+fluxCalFieldIds+"', # source name = "+fluxCalSourceName+"\n"
                                if fluxCalSourceName not in theFluxCalNames:
                                    theFluxCalNames.append(fluxCalSourceName)
                                casaCmd = casaCmd + "  spw = '"+str(spwIds[j])+"', # center frequency of spw = "+str(spwMeanFreq[j]/1.e9)+"GHz\n"
                                casaCmd = casaCmd + "  usescratch = True,\n" # because of bug in virt. model (see SCIREQ-2035)
                                casaCmd = casaCmd + "  standard = 'manual',\n"
                                casaCmd = casaCmd + "  fluxdensity = ["+str(flux1)+", 0, 0, 0]) # frequency of measurement = "+str(frequency1/1.e9)+"GHz\n\n"


    return casaCmd, list(theFluxCalNames)

#######################################

def doBandpassCalibration(msName, msName1='', bpassCalId='', chanAvg=1.0, refant='', iHaveSplitMyScienceSpw=False, 
                          calTableName=[], lowSNR=False, doplot=True, phaseDiff=False, solnorm=True, lbc=False, 
                          isB2B=None, isBWSW=None, combineB2BLFHFspws=None, isFullP=None, theFluxCalNames=[], vetoBPBootstrap=False, valueMaps={}):
    """Generate code for the bandpass calibration step of a calibration script."""

    print('\n*** doBandpassCalibration ***')
    print('Gathering information ...')

    casaCmd = ''

    if msName1 == '': 
        msName1 = msName
    if refant == '': 
        casalog.post('ERROR: No reference antenna specified.','SEVERE')
        return False
    if chanAvg > 1: 
        casalog.post('ERROR: The channel averaging bandwidth must be specified as a fraction of the total bandwidth.','SEVERE')
        return False
    if lowSNR == True: 
        chanAvg = 1.0

    if isB2B==None and isBWSW==None: # neither officially B2B nor BWSW but we double-check in case this method is called from outside the scriptgen
        isB2B, isBWSW = isB2BorBWSW(msName, valueMaps)

    if isB2B: 
        print('WARNING: TREATING THIS AS A B2B-TRANSFER OBSERVATION.')
        solnorm = True
        phaseDiff = True
        if combineB2BLFHFspws:
            print(' NOTE: WILL COMBINE RESPECTIVE LF AND HF SPWS')
    elif isBWSW:
        print('WARNING: TREATING THIS AS A BW-SWITCHING OBSERVATION.')
        print('WARNING: Forcing solnorm to False.')
        solnorm = False
        phaseDiff = True
    elif phaseDiff == True and solnorm == True:
        print('WARNING: phaseDiff True and not B2B, forcing solnorm to False.')
        solnorm = False

    if isFullP==None: # not officially full pol but we double-check in case this method is called from outside the scriptgen
        isFullP = isFullPol(msName, valueMaps)


    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/FIELD')
    fieldNames = mytb.getcol('NAME')
    mytb.close()

    ###

    print('Selecting bandpass calibrator ...')

    intentSources = sfsdr.getIntentsAndSourceNames(msName)

    if type(theFluxCalNames) != list or len(theFluxCalNames)==0:
        ampCalId = intentSources['CALIBRATE_AMPLI']['id'] + intentSources['CALIBRATE_FLUX']['id']
        ampCalId = np.unique([i for i in ampCalId if i != '']).tolist()

        if len(ampCalId) > 1:
            casaCmd = casaCmd + "# Note: there is more than one flux calibrator, picking the first one: "+fieldNames[ampCalId[0]]+".\n"
        ampCalId = ampCalId[0]
            
    else: # flux cal was given via theFluxCalNames
        ampCalId = list(fieldNames).index(theFluxCalNames[0])

    if bpassCalId == '':

        bpassCalId = intentSources['CALIBRATE_BANDPASS']['id']

        if bpassCalId[0] != '':

            if len(bpassCalId) != 1: casaCmd = casaCmd + "# Note: there is more than one bandpass calibrator, picking the first one: "+fieldNames[bpassCalId[0]]+".\n"
            bpassCalId = bpassCalId[0]

        elif isB2B:
            casaCmd = casaCmd + "# Note: there are no bandpass calibrators. Will try to use the diffgain calibrator instead.\n"
            phaseCalId = intentSources['CALIBRATE_DIFFGAIN']['id']
            if len(phaseCalId) > 0:
                bpassCalId = phaseCalId[0]
            else:
                casalog.post('ERROR: no bandpass or diffgain calibrator found in this alleged B2B observation.','SEVERE')
                return False
        else:
            casaCmd = casaCmd + "# Note: there are no bandpass calibrators. Will try to use a phase calibrator instead.\n"
            phaseCalId = intentSources['CALIBRATE_PHASE']['id']
            phaseOnlyCalId = [i for i in phaseCalId if i not in ampCalId]
            if len(phaseOnlyCalId) != 1: casaCmd = casaCmd + "# Note: there is more than one phase calibrator, picking the first one: "+fieldNames[phaseOnlyCalId[0]]+".\n"
            bpassCalId = phaseOnlyCalId[0]

    print('Will use field ID '+str(bpassCalId))

    bootstrap = False
    if bpassCalId != ampCalId:
        if not vetoBPBootstrap:
            casalog.post('The bandpass calibrator was not used as flux calibrator!', 'WARN')
            casalog.post('Will need to bootstrap its spectrum!', 'WARN')
            bootstrap = True
        else:
            casalog.post('The bandpass calibrator was not used as flux calibrator', 'WARN')
            casalog.post('but vetoBPBootstrap was set to True. So we will not bootstrap the BP spectrum!', 'WARN')
        
        
    ###

    spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
    spwIds = sorted(spwInfo.keys())

    sciNumChans = []
    sciChanWidths = []
    for i in spwIds: 
        sciNumChans.append(spwInfo[i]['numChans'])
        sciChanWidths.append(spwInfo[i]['meanAbsChanWidth'])
    if not isB2B and len(list(dict.fromkeys(sciNumChans).keys())) != 1:
        print('WARNING: This seems to be a mixed-mode dataset, please double-check the script.')

    spwSpec = ''
    for i in range(len(spwIds)):
        if spwSpec != '': spwSpec = spwSpec+','
        startChan = int((spwInfo[spwIds[i]]['numChans'] / 2.) * (1-chanAvg))
        endChan = int((spwInfo[spwIds[i]]['numChans'] / 2.) * (1+chanAvg))-1
        if iHaveSplitMyScienceSpw == True:
            spwSpec = spwSpec+str(i)
        else:
            spwSpec = spwSpec+str(spwIds[i])
        spwSpec = spwSpec+':'+str(startChan)+'~'+str(endChan)

    hasNoLFbpB2B = False
    if isB2B:
        spwInfoB = sfsdr.getSpwInfo(msName, intent='CALIBRATE_DIFFGAIN', caching=True)
        spwBIds = sorted(spwInfoB.keys())
        if spwBIds != spwIds:
            casalog.post('This B2B observation has no LF bandpass cal. Will use DiffGainCal.', 'WARN')
            hasNoLFbpB2B = True

            dgCalId = intentSources['CALIBRATE_DIFFGAIN']['id'][0]

            sciBNumChans = []
            for i in spwBIds: 
                sciBNumChans.append(spwInfoB[i]['numChans'])

            spwBSpec = ''
            for i in range(len(spwBIds)):
                if spwBSpec != '': spwBSpec = spwBSpec+','
                startChan = int((spwInfoB[spwBIds[i]]['numChans'] / 2.) * (1-chanAvg))
                endChan = int((spwInfoB[spwBIds[i]]['numChans'] / 2.) * (1+chanAvg))-1
                if iHaveSplitMyScienceSpw == True:
                    spwBSpec = spwBSpec+str(i)
                else:
                    spwBSpec = spwBSpec+str(spwBIds[i])
                spwBSpec = spwBSpec+':'+str(startChan)+'~'+str(endChan)

    ###

    print('Determining scan list ...')

    bpassCalScanList = []

    mymsmd = msmdtool()
    mymsmd.open(msName)
    
    if 'CALIBRATE_BANDPASS#ON_SOURCE' in mymsmd.intentsforfield(bpassCalId):
        bpassCalScanList = mymsmd.scansforintent('CALIBRATE_BANDPASS#ON_SOURCE') 
    else:
        bpassCalScanList = mymsmd.scansforfield(bpassCalId) 

    atmCalScanList = mymsmd.scansforintent('CALIBRATE_ATMOSPHERE*')
    bpassCalScanList = [i for i in bpassCalScanList if i not in atmCalScanList]  # still need index list for B2B filter - then str after B2B extra clause


    if isB2B:
        diffGainCalScanListLow  = ','.join([str(i) for i in mymsmd.scansforintent('CALIBRATE_DIFFGAIN#REFERENCE')])

        ### now we need additional things that are useful for B2B as the LF and HF are separated
        ### copied (and edited) code from Part I calibration. It is repeated also in part II.
        ### Could be rearranged if needed.


        ### determine the field, scan, and spw ids to be used in the code

        if msName in valueMaps.keys():
            vm = valueMaps[msName]
            print('Using canned ValueMap.')
        else:
            vm = aU.ValueMapping(msName)
            valueMaps[msName] = vm

        intentSources = sfsdr.getIntentsAndSourceNames(msName)
        diffGainCalId = intentSources['CALIBRATE_DIFFGAIN']['id'][0]

        spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#SIGNAL')
        if spwsDiffgainsig == []:
            spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#ON_SOURCE')
            sigintent = 'CALIBRATE_DIFFGAIN#ON_SOURCE'
        else:
            sigintent = 'CALIBRATE_DIFFGAIN#SIGNAL'


        spwHighDict = sfsdr.getSpwInfo(msName,intent=sigintent, caching=True)
        spwLowDict = sfsdr.getSpwInfo(msName,intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
        spwHigh = sorted(spwHighDict.keys())
        spwLow = sorted(spwLowDict.keys())


        # determine spw mapping
        if iHaveSplitMyScienceSpw == True:
            spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
            spwIds = sorted(spwInfo.keys())

            spwHighSplit = []
            for myspw in spwHigh:
                spwHighSplit.append(spwIds.index(myspw))  
            spwHighStr = ','.join([str(i) for i in spwHighSplit])

            spwLowSplit = []
            for myspw in spwLow:
                spwLowSplit.append(spwIds.index(myspw))
            spwLowStr = ','.join([str(i) for i in spwLowSplit])


            ## Mapping only needed LF to LF or HF to HF for Bandpass
            if combineB2BLFHFspws:
                LFtoLF_BP = list(range(max(spwLowSplit)+1))
                HFtoHF_BP = list(range(max(spwHighSplit)+1))
                for spwMap in spwLowSplit:
                    LFtoLF_BP[spwMap] = min(spwLowSplit)  # maps to lowest index
                for spwMap in spwHighSplit:
                    HFtoHF_BP[spwMap] = min(spwHighSplit)  # maps to highest index

                    

        else: # no reindexing took place
            spwHighStr = ','.join([str(i) for i in spwHigh])
            spwLowStr = ','.join([str(i) for i in spwLow])

            if combineB2BLFHFspws:
                LFtoLF_BP = list(range(max(spwLow)+1))
                HFtoHF_BP = list(range(max(spwHigh)+1))
                for spwMap in spwLow:
                    LFtoLF_BP[spwMap] = min(spwLow)  # maps to lowest index
                for spwMap in spwHigh:
                    HFtoHF_BP[spwMap] = min(spwHigh)  # maps to highest index


        # separate Bandpass scans
        # simply get from SpWs by exclusion
        
        bpassCalScanListLow =  [str(scnuse) for scnuse in bpassCalScanList if scnuse in mymsmd.scansforspw(min(spwLow))]
        bpassCalScanListHigh =  [str(scnuse) for scnuse in  bpassCalScanList if scnuse in mymsmd.scansforspw(min(spwHigh))]
        bpassCalScanListLow = ','.join(bpassCalScanListLow)   
        bpassCalScanListHigh = ','.join(bpassCalScanListHigh)  

        
    ### end B2B extras  ###
        
    mymsmd.close()

    ### now make the string list and print (global script gen) - filtered ATM out above isB2B 
    
    bpassCalScanList = [str(i) for i in bpassCalScanList]

    print(bpassCalScanList)

    bpassCalScanList = ','.join(bpassCalScanList)  

    #######
    if isB2B and combineB2BLFHFspws:
        print('bpassCalScanListLow ', bpassCalScanListLow)
        print('bpassCalScanListHigh ', bpassCalScanListHigh)
        print('spwHighStr ', spwHighStr)
        print('spwLowStr ', spwLowStr)
        print('LFtoLF_BP ', str(LFtoLF_BP))
        print('HFtoHF_BP ', str(HFtoHF_BP))

    ######

    print('Writing code ...')

    calTableName1 = msName1+'.bandpass'

    if isB2B:
       # permanent phasediff - but not in same sense as implemented below (could be useful for all modes)
        casaCmd = casaCmd + "os.system('rm -rf %s.ap_pre_phasediff') \n"%(msName1)
        casaCmd = casaCmd + "\ngaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".ap_pre_phasediff',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
        casaCmd = casaCmd + "  spw = '"+spwSpec+"',\n"
        if bpassCalScanList != '': 
            casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"
        casaCmd = casaCmd + "  solint = 'inf',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        casaCmd = casaCmd + "  refantmode = 'strict',\n"
        casaCmd = casaCmd + "  minsnr = 3.0,\n"
        casaCmd = casaCmd + "  calmode = 'p')\n"

        if hasNoLFbpB2B:  ## likely defunct, only useful for old non-conforming data
            casaCmd = casaCmd + "\ngaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+msName1+".ap_pre_phasediff',\n"
            casaCmd = casaCmd + "  field = '"+str(dgCalId)+"', # "+fieldNames[int(dgCalId)]+"\n"
            casaCmd = casaCmd + "  spw = '"+spwBSpec+"',\n"
            casaCmd = casaCmd + "  scan = '"+diffGainCalScanListLow+"',\n"
            casaCmd = casaCmd + "  solint = 'inf',\n"
            casaCmd = casaCmd + "  minsnr = 3.0,\n"
            casaCmd = casaCmd + "  refant = '"+refant+"',\n"
            casaCmd = casaCmd + "  refantmode = 'strict',\n"
            casaCmd = casaCmd + "  calmode = 'p',\n"
            casaCmd = casaCmd + "  append = True)\n"


        ##  now we do the ap pre bandpass, and must always use the ap_pre_phasediff in solve    
            
        casaCmd = casaCmd + "\nos.system('rm -rf %s.ap_pre_bandpass') \n"%(msName1)
        casaCmd = casaCmd + "\ngaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".ap_pre_bandpass',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
        casaCmd = casaCmd + "  spw = '"+spwSpec+"',\n"
        if bpassCalScanList != '': 
            casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"
        casaCmd = casaCmd + "  solint = 'int',\n"
        if combineB2BLFHFspws:
            casaCmd = casaCmd + "  combine='spw',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        casaCmd = casaCmd + "  refantmode = 'strict',\n"
        casaCmd = casaCmd + "  minsnr = 3.0,\n"
        casaCmd = casaCmd + "  gaintable = ['%s.ap_pre_phasediff'], \n"%(msName1)
        casaCmd = casaCmd + "  calmode = 'p')\n"

        if hasNoLFbpB2B:
            casaCmd = casaCmd + "\ngaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+msName1+".ap_pre_bandpass',\n"
            casaCmd = casaCmd + "  field = '"+str(dgCalId)+"', # "+fieldNames[int(dgCalId)]+"\n"
            casaCmd = casaCmd + "  spw = '"+spwBSpec+"',\n"
            casaCmd = casaCmd + "  scan = '"+diffGainCalScanListLow+"',\n"
            casaCmd = casaCmd + "  solint = 'int',\n"
            casaCmd = casaCmd + "  minsnr = 3.0,\n"
            casaCmd = casaCmd + "  refant = '"+refant+"',\n"
            casaCmd = casaCmd + "  refantmode = 'strict',\n"
            casaCmd = casaCmd + "  calmode = 'p',\n"
            casaCmd = casaCmd + "  gaintable = ['%s.ap_pre_phasediff'], \n"%(msName1)
            casaCmd = casaCmd + "  append = True)\n"

 

    else:  ## i.e. NOT B2B, as it was before
        casaCmd = casaCmd + "os.system('rm -rf %s.ap_pre_bandpass') \n"%(msName1)
        casaCmd = casaCmd + "\ngaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".ap_pre_bandpass',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
        casaCmd = casaCmd + "  spw = '"+spwSpec+"',\n"
        if bpassCalScanList != '': 
            casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"
        casaCmd = casaCmd + "  solint = 'int',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        if isFullP:   # removed B2B option as it is a stand alone above
            casaCmd = casaCmd + "  refantmode = 'strict',\n"
        casaCmd = casaCmd + "  calmode = 'p')\n"


    if doplot:
        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".ap_pre_bandpass', msName='"+msName1+"', interactive=False) \n\n"

    if bootstrap and not phaseDiff: # bootstrap the bandpass calibrator spectrum (SCIREQ-2477)

        casaCmd = casaCmd + "os.system('rm -rf "+msName1+".bandpass_tentative')\n" 
        casaCmd = casaCmd + "bandpass(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".bandpass_tentative',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
        if bpassCalScanList != '': 
            casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"
            
        casaCmd = casaCmd + "  solint = 'inf',\n"
        casaCmd = casaCmd + "  combine = 'scan',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        casaCmd = casaCmd + "  solnorm = "+str(solnorm)+",\n"
        casaCmd = casaCmd + "  bandtype = 'B',\n"
        casaCmd = casaCmd + "  gaintable = '"+msName1+".ap_pre_bandpass')\n\n"

        if doplot:
            casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".bandpass_tentative', msName='"+msName1+"', interactive=False)\n\n" 
    
        casaCmd = casaCmd + "# bootstrapping BP calibrator spectrum\n\n"

        casaCmd = casaCmd + "os.system('rm -rf "+msName1+".G1p')\n"
        casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".G1p',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+","+str(ampCalId)+"', # bandpass, flux calibrator\n"
        casaCmd = casaCmd + "  gaintable = '"+msName1+".bandpass_tentative',\n"
        casaCmd = casaCmd + "  gaintype = 'G',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        casaCmd = casaCmd + "  calmode = 'p',\n"
        casaCmd = casaCmd + "  solint = 'int')\n\n"

        casaCmd = casaCmd + "os.system('rm -rf "+msName1+".G1')\n"
        casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".G1',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+","+str(ampCalId)+"', # bandpass, flux calibrator\n"
        casaCmd = casaCmd + "  gaintable = ['"+msName1+".bandpass_tentative', '"+msName1+".G1p'],\n"
        casaCmd = casaCmd + "  gaintype = 'G',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        casaCmd = casaCmd + "  calmode = 'ap',\n"
        casaCmd = casaCmd + "  solint = 'inf')\n\n"

        casaCmd = casaCmd + "os.system('rm -rf "+msName1+".F1')\n"
        casaCmd = casaCmd + "flux_bandpass = fluxscale(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".G1',\n"
        casaCmd = casaCmd + "  fluxtable = '"+msName1+".F1',\n"
        casaCmd = casaCmd + "  reference = '"+str(ampCalId)+"',\n"
        casaCmd = casaCmd + "  transfer = '"+str(bpassCalId)+"',\n"
        casaCmd = casaCmd + "  listfile = '"+msName1+".bandpass.fluxinfo',\n"
        casaCmd = casaCmd + "  fitorder=1)\n\n"
        
        casaCmd = casaCmd + "setjy(vis='"+msName1+"',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"',\n"
        casaCmd = casaCmd + "  scalebychan = True,\n"
        casaCmd = casaCmd + "  standard = 'fluxscale',\n"
        casaCmd = casaCmd + "  fluxdict = flux_bandpass)\n\n"
        
        casaCmd = casaCmd + "os.system('rm -rf "+msName1+".G0.b')\n"
        casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".G0.b',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # bandpass\n"
        casaCmd = casaCmd + "  gaintype = 'G',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        casaCmd = casaCmd + "  calmode = 'p',\n"
        casaCmd = casaCmd + "  solint = 'int')\n\n"

        casaCmd = casaCmd + "os.system('rm -rf "+calTableName1+"')\n"
        casaCmd = casaCmd + "bandpass(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+calTableName1+"',\n"
        casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # bandpass\n"
        if bpassCalScanList != '': 
            casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"

        casaCmd = casaCmd + "  solint = 'inf',\n"
        casaCmd = casaCmd + "  combine = 'scan',\n"
        casaCmd = casaCmd + "  refant = '"+refant+"',\n"
        casaCmd = casaCmd + "  solnorm = True,\n"
        casaCmd = casaCmd + "  bandtype = 'B',\n"
        casaCmd = casaCmd + "  gaintable = '"+msName1+".G0.b')\n\n"

        if doplot:
            casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+"', msName='"+msName1+"', interactive=False) \n"

        calTableName.append(calTableName1)
        
    else: # no bootstrapping   

        if bootstrap: # bootstrapping was not done because phasediff was True at the same time
            casalog.post('WARNING: bandpass cal != fluxcal but phaseDiff mode is on. No bandpass bootstrapping implemented.', 'WARN')

        ## FOR B2B only LF and HF separated 
        if isB2B:
            # doing the LF first
            casaCmd = casaCmd + "os.system('rm -rf %s.bandpass') \n"%(msName1)

            if not hasNoLFbpB2B:
                casaCmd = casaCmd + "bandpass(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+calTableName1+"',\n"
                casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
                if bpassCalScanList != '': 
                    casaCmd = casaCmd + "  scan = '"+bpassCalScanListLow+"',\n"  
                casaCmd = casaCmd + "  solint = 'inf',\n"
                casaCmd = casaCmd + "  minsnr = 3.0,\n"
                casaCmd = casaCmd + "  combine = 'scan',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                casaCmd = casaCmd + "  solnorm = "+str(solnorm)+",\n"
                casaCmd = casaCmd + "  bandtype = 'B',\n"
                if combineB2BLFHFspws:  ## ADD THE MAPS
                    casaCmd = casaCmd + "  spwmap = [[],"+str(LFtoLF_BP)+"], \n"  
                casaCmd = casaCmd + "  gaintable = ['"+msName1+".ap_pre_phasediff','"+msName1+".ap_pre_bandpass'])\n"
 
            #if hasNoLFbpB2B:
            else:
                casaCmd = casaCmd + "bandpass(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+calTableName1+"',\n"
                casaCmd = casaCmd + "  field = '"+str(dgCalId)+"', # "+fieldNames[int(dgCalId)]+"\n"
                casaCmd = casaCmd + "  scan = '"+diffGainCalScanListLow+"',\n"
                casaCmd = casaCmd + "  solint = 'inf',\n"
                casaCmd = casaCmd + "  minsnr = 3.0,\n"
                casaCmd = casaCmd + "  combine = 'scan',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                casaCmd = casaCmd + "  solnorm = "+str(solnorm)+",\n"
                casaCmd = casaCmd + "  bandtype = 'B',\n"
                #casaCmd = casaCmd + "  append = True,\n"
                if combineB2BLFHFspws:  ## ADD THE MAPS - don't know if this functions given no LF bandpass and other field used - old data, wont use? 
                    casaCmd = casaCmd + "  spwmap = [[],"+str(LFtoLF_BP)+"], \n"  
                casaCmd = casaCmd + "  gaintable = ['"+msName1+".ap_pre_phasediff','"+msName1+".ap_pre_bandpass'])\n"

            ## here wanted for the HF wanted a run over AU tool
            ## but that needs the MS and the Tsys table
            ## cannnot be calling aU in the script as users do not have it by default
            ## placeholder for now just to add instructins to analysts

            casaCmd = casaCmd + "\n  # Note for analyst, the below needs to be copied/edited \n"
            casaCmd = casaCmd + "  # possibly for each SpW individaully to apply X channel binning as required \n"
            casaCmd = casaCmd + "  # PLEASE delete this comment before delivering the script \n"            
            casaCmd = casaCmd + "bandpass(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+calTableName1+"',\n"
            casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
            if bpassCalScanList != '': 
                casaCmd = casaCmd + "  scan = '"+bpassCalScanListHigh+"',\n"
            casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n" 
            casaCmd = casaCmd + "  solint = 'inf,Xch',\n"
            casaCmd = casaCmd + "  minsnr = 3.0,\n"
            casaCmd = casaCmd + "  combine = 'scan',\n"
            casaCmd = casaCmd + "  refant = '"+refant+"',\n"
            casaCmd = casaCmd + "  solnorm = "+str(solnorm)+",\n"
            casaCmd = casaCmd + "  bandtype = 'B',\n"
            casaCmd = casaCmd + "  append = True,\n"
            if combineB2BLFHFspws:  ## ADD THE MAPS
                    casaCmd = casaCmd + "  spwmap = [[],"+str(HFtoHF_BP)+"], \n"  
            casaCmd = casaCmd + "  gaintable = ['"+msName1+".ap_pre_phasediff','"+msName1+".ap_pre_bandpass'])\n"         
                
            if doplot:
                casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+"', msName='"+msName1+"', interactive=False) \n"


        else: #STANDARD          
            casaCmd = casaCmd + "os.system('rm -rf %s.bandpass') \n"%(msName1)
            casaCmd = casaCmd + "bandpass(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+calTableName1+"',\n"
            casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
            if bpassCalScanList != '': 
                casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"  
            casaCmd = casaCmd + "  solint = 'inf',\n"
            casaCmd = casaCmd + "  combine = 'scan',\n"
            casaCmd = casaCmd + "  refant = '"+refant+"',\n"
            casaCmd = casaCmd + "  solnorm = "+str(solnorm)+",\n"
            casaCmd = casaCmd + "  bandtype = 'B',\n"
            casaCmd = casaCmd + "  gaintable = '"+msName1+".ap_pre_bandpass')\n"

            if doplot:
                casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+"', msName='"+msName1+"', interactive=False) \n"

            minSciNumChans = min(sciNumChans)
            minSciChanWidth = min(sciChanWidths)

            if (minSciNumChans > 256 or isBWSW) and minSciChanWidth < 8E6:
                casaCmd = casaCmd + "\nos.system('rm -rf %s.bandpass_smooth20ch') \n"%(msName1)
                casaCmd = casaCmd + "\nbandpass(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+calTableName1+'_smooth20ch'+"',\n"
                casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
                if bpassCalScanList != '': 
                    casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"
                if lbc or isBWSW:
                    casaCmd = casaCmd + "  solint = 'inf,8MHz',\n"
                else:
                    casaCmd = casaCmd + "  solint = 'inf,20ch',\n"
                casaCmd = casaCmd + "  combine = 'scan',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                casaCmd = casaCmd + "  solnorm = "+str(solnorm)+",\n"
                casaCmd = casaCmd + "  bandtype = 'B',\n"
                casaCmd = casaCmd + "  gaintable = '"+msName1+".ap_pre_bandpass')\n"


                if doplot:
                    casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+'_smooth20ch'+"', msName='"+msName1+"', interactive=False) \n"

                #use the smoothed table from now on
                calTableName1 = calTableName1+'_smooth20ch'


        calTableName.append(calTableName1)

        if phaseDiff == True:
            if bpassCalId != ampCalId:

                casaCmd = casaCmd + "\n\n"

                fluxscaleDictName = []
                casaCmd = casaCmd + doGainCalibration(msName, msName1=msName1, refant=refant, bandpass=calTableName1, calmode2='a', phaseDiff=False, fluxscaleDictName=fluxscaleDictName, iHaveSplitMyScienceSpw=iHaveSplitMyScienceSpw, isBWSW=False, isFullP=isFullP, valueMaps=valueMaps)

                casaCmd = casaCmd + "\nsetjy(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
                casaCmd = casaCmd + "  standard = 'manual',\n"
                casaCmd = casaCmd + "  spw = '',\n"
                casaCmd = casaCmd + "  fluxdensity = "+fluxscaleDictName[0]+"['"+str(bpassCalId)+"']['fitFluxd'],\n"
                casaCmd = casaCmd + "  spix = "+fluxscaleDictName[0]+"['"+str(bpassCalId)+"']['spidx'][1],\n"  # Added trailing [1] - T. Hunter 2014-08-11
                casaCmd = casaCmd + "  reffreq = '%fGHz'%(1e-9*"+fluxscaleDictName[0]+"['"+str(bpassCalId)+"']['fitRefFreq']))\n\n" # T. Hunter 2014-08-11

                calTableName1 = msName1+'.bandpass2'
                casaCmd = casaCmd + "os.system('rm -rf %s.bandpass2') \n"%(msName1) # Added by CLB
                casaCmd = casaCmd + "bandpass(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+calTableName1+"',\n"
                casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
                if bpassCalScanList != '': 
                    casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"
                casaCmd = casaCmd + "  solint = 'inf',\n"
                casaCmd = casaCmd + "  combine = 'scan',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                casaCmd = casaCmd + "  solnorm = False,\n"
                casaCmd = casaCmd + "  bandtype = 'B',\n"
                casaCmd = casaCmd + "  gaintable = '"+msName1+".ap_pre_bandpass')\n"

                if minSciNumChans > 256 and minSciChanWidth < 8E6:
                    casaCmd = casaCmd + "\nos.system('rm -rf %s.bandpass2_smooth20ch') \n"%(msName1)
                    casaCmd = casaCmd + "\nbandpass(vis = '"+msName1+"',\n"
                    casaCmd = casaCmd + "  caltable = '"+calTableName1+'_smooth20ch'+"',\n"
                    casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[int(bpassCalId)]+"\n"
                    if bpassCalScanList != '': 
                        casaCmd = casaCmd + "  scan = '"+bpassCalScanList+"',\n"
                    casaCmd = casaCmd + "  solint = 'inf,20ch',\n"
                    casaCmd = casaCmd + "  combine = 'scan',\n"
                    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                    casaCmd = casaCmd + "  solnorm = False,\n"
                    casaCmd = casaCmd + "  bandtype = 'B',\n"
                    casaCmd = casaCmd + "  gaintable = '"+msName1+".ap_pre_bandpass')\n"
                    if doplot == True:
                        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+'_smooth20ch'+"', msName='"+msName1+"', interactive=False) \n"

                if doplot == True:
                    casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+"', msName='"+msName1+"', interactive=False) \n"

                if minSciNumChans > 256 and minSciChanWidth < 8E6:
                    calTableName1 = calTableName1+'_smooth20ch'

                calTableName[0] = calTableName1
                
        # endif phasediff        
    # endif bootstrap
                
    return casaCmd

##################################

def doGainCalibration(msName, msName1='', refant='', bandpass='', gaintypeForAmp='T', doplot=True, calFieldsOnly=True, calmode2='ap', phaseDiff='', phaseDiffCalTableName=[], fluxscaleDictName=[], ampForSci=[], iHaveSplitMyScienceSpw=False, isBWSW=None, isFullP=None, valueMaps={}):
    """Generate code for the gain calibration step of a calibration script."""

    print('\n*** doGainCalibration ***')
    print('Gathering information ...')

    if msName1 == '': msName1 = msName
    if refant == '': 
        casalog.post('ERROR: No reference antenna specified.', 'SEVERE')
        return False
    if bandpass == '': 
        casalog.post('ERROR: No bandpass cal table specified.', 'SEVERE')
        return False

    myCasaVersion = aU.getCasaVersion()

    fluxscaleDictName.append('fluxscaleDict')

    print('Selecting phase calibrator(s) ...')

    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/FIELD')
    fieldNames = mytb.getcol('NAME')
    mytb.close()

    fieldIds = list(range(len(fieldNames)))

    intentSources = sfsdr.getIntentsAndSourceNames(msName)
    sciFieldIds = intentSources['OBSERVE_TARGET']['id']
    phaseCalIds = intentSources['CALIBRATE_PHASE']['id']
    print('sciFieldIds', sciFieldIds)
    print('phaseCalIds', phaseCalIds)

    if 'OBSERVE_CHECK' in list(intentSources.keys()): 
        sciFieldIds += intentSources['OBSERVE_CHECK']['id']
    if sciFieldIds[0] == '': 
        casalog.post('THERE SEEMS TO BE NO OBSERVE_TARGET OR OBSERVE_CHECK FIELD', 'WARN')
    if phaseCalIds[0] == '': 
        casalog.post('THERE SEEMS TO BE NO CALIBRATE_PHASE FIELD', 'WARN')
    else:
        for i in phaseCalIds:
            if i in sciFieldIds:
                casalog.post('THE CALIBRATE_PHASE FIELD '+str(i)+' ALSO HAS A TARGET OR CHECK INTENT!', 'WARN')
                casalog.post('IT CAN THEREFORE BY DEFAULT NOT BE USED AS A PHASE CALIBRATOR!', 'WARN')
                casalog.post('Please investigate and rectify the calibration script draft by hand if needed.', 'WARN')

    if calFieldsOnly:
        calFieldIds = [i for i in fieldIds if i not in sciFieldIds]
    else:
        calFieldIds = fieldIds

    mymsmd = msmdtool()
    mymsmd.open(msName)
    hasdata = []
    for i in calFieldIds:
        calFieldIntents = mymsmd.intentsforfield(i)
        hasdata1 = 0
        for j in calFieldIntents:
            if re.search('^CALIBRATE_(POINTING|ATMOSPHERE|WVR)', j) == None:
                hasdata1 = 1
                break
        hasdata.append(hasdata1)
    mymsmd.close()

    calFieldIds = [calFieldIds[i] for i in range(len(calFieldIds)) if hasdata[i] == 1]

    print('calFieldIds', calFieldIds)

    if len(calFieldIds) == '': 
        casalog.post('ERROR: There seems to be no calibrator field.', 'SEVERE')
        return False

    calFieldNames = [fieldNames[i] for i in calFieldIds]
    calFieldNames = ','.join(calFieldNames)

    if len(calFieldIds) > 1:
        j0 = 0
        calFieldIds1 = str(calFieldIds[j0])
        for j in range(len(calFieldIds)-1):
            if calFieldIds[j+1] == calFieldIds[j]+1: continue
            calFieldIds1 = calFieldIds1 + '~' + str(calFieldIds[j])
            j0 = j+1
            calFieldIds1 = calFieldIds1 + ',' + str(calFieldIds[j0])
        calFieldIds1 = calFieldIds1 + '~' + str(calFieldIds[j+1])
        print('Multiple phase calibrators will be used: ', calFieldIds1)
    else:
        calFieldIds1 = str(calFieldIds[0])
        print('Single phase calibrator: ', calFieldIds1)

    spwInfo = sfsdr.getSpwInfo(msName, caching=True)
    spwIds = sorted(spwInfo.keys())

    ###

    if isBWSW==None: # not officially BWSW but we double-check in case this method is called from outside the scriptgen
        isB2B, isBWSW = isB2BorBWSW(msName, valueMaps)
        if isB2B: 
            casalog.post('This is a B2B observation. doGainCalibration does not support that. Use doB2BGainCalibrationPartI() and ...PartII() instead of doGainCalibration()', 'SEVERE')
            return False

    if isBWSW:
        casalog.post('TREATING THIS AS A BW-SWITCHING OBSERVATION.','WARN')
        print('WARNING: Forcing phaseDiff=True.')
        phaseDiff = True

    if isFullP==None: # not officially full pol but we double-check in case this method is called from outside the scriptgen
        isFullP = isFullPol(msName, valueMaps)

    ###

    casaCmd = ''

    if phaseDiff == True:

        if sciFieldIds[0] != '':

            bpassCalId = intentSources['CALIBRATE_BANDPASS']['id']
            if bpassCalId[0] == '': 
                casalog.post('ERROR: There is no bandpass calibrator.', 'SEVERE')
                return False
            if len(bpassCalId) != 1: 
                casaCmd = casaCmd + "# Note: there is more than one bandpass calibrator, I'm picking the first one: "+fieldNames[bpassCalId[0]]+".\n"
            bpassCalId = bpassCalId[0]

            diffGainCalId = ''
            diffGainCalScanList = ''
            if 'CALIBRATE_DIFFGAIN' in list(intentSources.keys()):
                diffGainCalId = intentSources['CALIBRATE_DIFFGAIN']['id']
                if len(diffGainCalId) != 1: 
                    casaCmd = casaCmd + "# Note: there is more than one diffgain calibrator, I'm picking the first one: "+fieldNames[diffGainCalId[0]]+".\n"
                diffGainCalId = diffGainCalId[0]
            if diffGainCalId == '':
                diffGainCalId = bpassCalId
            else:
                mymsmd = msmdtool()
                mymsmd.open(msName)
                diffGainCalScanList = mymsmd.scansforintent('CALIBRATE_DIFFGAIN#*').tolist()
                mymsmd.close()
                diffGainCalScanList = [str(j) for j in diffGainCalScanList]
                diffGainCalScanList = ','.join(diffGainCalScanList)

            casaCmd = casaCmd + "os.system('rm -rf %s.phasediff_inf') \n"%(msName1)
            casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+msName1+".phasediff_inf',\n"
            casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
            if diffGainCalScanList != '': 
                casaCmd = casaCmd + "  scan = '"+diffGainCalScanList+"',\n"
            casaCmd = casaCmd + "  solint = 'inf',\n"
            casaCmd = casaCmd + "  combine = 'scan',\n"
            casaCmd = casaCmd + "  refant = '"+refant+"',\n"
            if isFullP:
                casaCmd = casaCmd + "  refantmode = 'strict',\n"
            casaCmd = casaCmd + "  gaintype = 'G',\n"
            casaCmd = casaCmd + "  calmode = 'p',\n"
            casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n\n"

            phaseDiffCalTableName.append(msName1+'.phasediff_inf')

            if doplot == True: 
                casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".phasediff_inf', msName='"+msName1+"', interactive=False) \n\n"

            casaCmd = casaCmd + "for i in "+str(calFieldIds)+": # "+calFieldNames+"\n"
            casaCmd = casaCmd + "  os.system('rm -rf %s.phase_int'+str(i)) \n"%(msName1)
            casaCmd = casaCmd + "  gaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "    caltable = '"+msName1+".phase_int'+str(i),\n"
            casaCmd = casaCmd + "    field = str(i),\n"
            casaCmd = casaCmd + "    solint = 'int',\n"
            casaCmd = casaCmd + "    #combine = 'spw', # change calspwmap below and also combine for phase_inf if you uncomment this\n"
            casaCmd = casaCmd + "    refant = '"+refant+"',\n"
            if isFullP:
                casaCmd = casaCmd + "    refantmode = 'strict',\n"
            casaCmd = casaCmd + "    gaintype = 'G',\n"
            casaCmd = casaCmd + "    calmode = 'p',\n"
            casaCmd = casaCmd + "    append = False,\n"
            casaCmd = casaCmd + "    gaintable = ['"+bandpass+"', '"+msName1+".phasediff_inf'])\n\n"

            if doplot == True: 
                casaCmd = casaCmd + "  if applyonly != True: es.checkCalTable('"+msName1+".phase_int'+str(i), msName='"+msName1+"', interactive=False) \n\n"

            ###
            # calspwmap creation for spw-combined and non-spw-combined phase_int

            mymsmd = msmdtool()
            mymsmd.open(msName)

            spwInfo3 = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
            spwIds3 = sorted(spwInfo3.keys())
            spwInfoDGCR = sfsdr.getSpwInfo(msName, intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
            spwIdsDGCR = sorted(spwInfoDGCR.keys())
            spwInfoScience = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET#ON_SOURCE', caching=True)
            spwIdsScience = sorted(spwInfoScience.keys())

            if iHaveSplitMyScienceSpw == True:
                numSpws = len(spwIds3)
                tmpIds = []
                for i in spwIdsDGCR:
                    tmpIds.append(spwIds3.index(i))
                spwIdsDGCR = sorted(tmpIds)
                tmpIds = []
                for i in spwIdsScience:
                    tmpIds.append(spwIds3.index(i))
                spwIdsScience = sorted(tmpIds)
                tmpIds = []
                for i in spwIds3:
                    tmpIds.append(spwIds3.index(i))
                spwIds3 = sorted(tmpIds)

            else:
                numSpws = max(spwIds3)+2

            calspwmapcomb = [] # version to be used with spw-combination

            # decide whether the narrow or the wide SPWs come first
            narrowFirst = True
            for i in spwIdsScience:
                for j in spwIdsDGCR: 
                    if i > j:
                        narrowFirst = False
                        break

            # create calspwmapcomb accordingly
            minDGCRSpw = min(spwIdsDGCR)
            minSciSpw = min(spwIdsScience)
            if narrowFirst:
                for i in range(minDGCRSpw):
                    calspwmapcomb.append(minSciSpw)
                for i in range(numSpws-minDGCRSpw):
                    calspwmapcomb.append(minDGCRSpw)
            else:
                for i in range(minSciSpw):
                    calspwmapcomb.append(minDGCRSpw)
                for i in range(numSpws-minSciSpw):
                    calspwmapcomb.append(minSciSpw)


            casaCmd = casaCmd + "# Case combine='spw' in phase_int\n"
            casaCmd = casaCmd + "#calspwmap = "+repr(calspwmapcomb)+"\n\n"
            casaCmd = casaCmd + "# Case combine='' in phase_int\n"
            casaCmd = casaCmd + "calspwmap = list(range("+str(numSpws)+"))\n\n"

            ## create calspwmap for flux scale
            calspwmapf = []

            for i in range(numSpws):
                if i in spwIdsScience:
                    myDGCRSpw = spwIdsDGCR[mymsmd.baseband(i)-1] # wide SPWs are ordered by baseband number
                    calspwmapf.append(myDGCRSpw)
                else:
                    calspwmapf.append(i)

            mymsmd.close()        

            ###

            casaCmd = casaCmd + "os.system('rm -rf %s.ampli_inf') \n"%(msName1)
            casaCmd = casaCmd + "for i in "+str(calFieldIds)+": # "+calFieldNames+"\n"
            if myCasaVersion > '6.3.0': # accommodate changed gaincal append behaviour in 6.4
                casaCmd = casaCmd + "  myappend = False if i == "+str(calFieldIds[0])+" else True\n"
            casaCmd = casaCmd + "  gaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "    caltable = '"+msName1+".ampli_inf',\n"
            casaCmd = casaCmd + "    field = str(i),\n"
            casaCmd = casaCmd + "    solint = 'inf',\n"
            casaCmd = casaCmd + "    refant = '"+refant+"',\n"
            if isFullP:
                casaCmd = casaCmd + "    refantmode = 'strict',\n"
            casaCmd = casaCmd + "    gaintype = '"+gaintypeForAmp+"',\n"
            casaCmd = casaCmd + "    calmode = 'a',\n"
            if myCasaVersion > '6.3.0': # accommodate changed gaincal append behaviour in 6.4
                casaCmd = casaCmd + "    append = myappend,\n"
            else:
                casaCmd = casaCmd + "    append = True,\n"
            casaCmd = casaCmd + "    gaintable = ['"+bandpass+"', '"+msName1+".phasediff_inf', '"+msName1+".phase_int'+str(i)],\n"
            casaCmd = casaCmd + "    spwmap = [[], [], calspwmap])\n\n"

            ampForSci.append(msName1+'.ampli_inf')

            if doplot == True: 
                casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".ampli_inf', msName='"+msName1+"', interactive=False) \n\n"

            casaCmd = casaCmd + "# spwmap for flux inf: matched by baseband to preserve signal path and keep science SPWs statistically independent as in non-BWSW case\n"
            casaCmd = casaCmd + "calspwmapf = "+repr(calspwmapf)+"\n\n"

            casaCmd = casaCmd + "os.system('rm -rf %s.flux_inf') \n"%(msName1)
            casaCmd = casaCmd + "os.system('rm -rf %s.fluxscale') \n"%(msName1)
            casaCmd = casaCmd + "mylogfile = casalog.logfile()\n"
            casaCmd = casaCmd + "casalog.setlogfile('"+msName1+".fluxscale')\n\n"
            casaCmd = casaCmd + fluxscaleDictName[0] + " = fluxscale(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"
            casaCmd = casaCmd + "  fluxtable = '"+msName1+".flux_inf',\n"
            casaCmd = casaCmd + "  reference = '"+str(bpassCalId)+"', # "+fieldNames[bpassCalId]+"\n"
            casaCmd = casaCmd + "  refspwmap = calspwmapf,\n"
            casaCmd = casaCmd + "  incremental = True)\n\n"
            casaCmd = casaCmd + "casalog.setlogfile(mylogfile)\n\n"
            casaCmd = casaCmd + "if applyonly != True: es.fluxscale2(caltable = '"+msName1+".ampli_inf', removeOutliers=True, msName='"+msName+"', writeToFile=True, preavg=10000)\n\n"

            casaCmd = casaCmd + "for i in "+str(phaseCalIds)+": # phasecal\n"
            casaCmd = casaCmd + "  os.system('rm -rf "+msName1+".phase_inf'+str(i)) \n"
            casaCmd = casaCmd + "  gaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "    caltable = '"+msName1+".phase_inf'+str(i),\n"
            casaCmd = casaCmd + "    field = str(i),\n"
            casaCmd = casaCmd + "    solint = 'inf',\n"
            casaCmd = casaCmd + "    #combine = 'spw',\n"
            casaCmd = casaCmd + "    refant = '"+refant+"',\n"
            if isFullP:
                casaCmd = casaCmd + "    refantmode = 'strict',\n"
            casaCmd = casaCmd + "    gaintype = 'G',\n"
            casaCmd = casaCmd + "    calmode = 'p',\n"
            casaCmd = casaCmd + "    append = False,\n"
            casaCmd = casaCmd + "    gaintable = ['"+bandpass+"', '"+msName1+".phasediff_inf'])\n\n"

            if doplot == True: 
                casaCmd = casaCmd + "  if applyonly != True: es.checkCalTable('"+msName1+".phase_inf'+str(i), msName='"+msName1+"', interactive=False) \n"
        # endif there are scifields

    else: # phaseDiff != True
        fluxCalId = sfsdr.getFieldsForSetjy(msName)

        if len(fluxCalId) == 0:
            casaCmd = casaCmd + "os.system('rm -rf %s.phase_int') \n"%(msName1)
            casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+msName1+".phase_int',\n"
            casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
            casaCmd = casaCmd + "  solint = 'int',\n"
            casaCmd = casaCmd + "  refant = '"+refant+"',\n"
            if isFullP:
                casaCmd = casaCmd + "  refantmode = 'strict',\n"
            casaCmd = casaCmd + "  gaintype = 'G',\n"
            casaCmd = casaCmd + "  calmode = 'p',\n"
            casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n\n"

            if doplot == True: 
                casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".phase_int', msName='"+msName1+"', interactive=False) \n\n"

            intentSources1 = intentSources['CALIBRATE_FLUX']

            if len(intentSources1['sourceid']) != 0:

                sourceFluxes = sfsdr.getFluxesFromSourceTable(msName)

                fluxCalSourceId = intentSources1['sourceid']

                if len(fluxCalSourceId) > 1:
                    casalog.post("THERE IS MORE THAN ONE FLUX CALIBRATOR. I WILL PICK THE FIRST ONE. THIS MAY BE WRONG.", 'WARN')
                    fluxCalSourceId = [j for j in fluxCalSourceId if j in list(sourceFluxes.keys())]

                if len(fluxCalSourceId) == 0: 
                    casalog.post("ERROR: There is no flux calibrator.", 'SEVERE')
                    return False

                fluxCalSourceId = fluxCalSourceId[0]

                fluxCalSourceName = intentSources1['name'][intentSources1['sourceid'].index(fluxCalSourceId)]

                if len(intentSources1['sourceid']) > 1:
                    mytb.open(msName+'/FIELD')
                    tb1 = mytb.query('SOURCE_ID == '+str(fluxCalSourceId))
                    fluxCalFieldIds1 = tb1.rownumbers().tolist()
                    tb1.close()
                    mytb.close()
                    fluxCalFieldIds = [j for j in fluxCalFieldIds1 if j in intentSources1['id']]
                else:
                    fluxCalFieldIds = intentSources1['id']

                if fluxCalSourceId in sourceFluxes:

                    if fluxCalSourceName != sourceFluxes[fluxCalSourceId]['sourceName']: 
                        casalog.post("ERROR: Source names do not match.", 'SEVERE')
                        return False

                    fluxCalId = fluxCalFieldIds[:]

            if len(fluxCalId) != 0:

                if len(fluxCalId) > 1: 
                    casaCmd = casaCmd + "# Note: there is more than one flux calibrator in this dataset, I'm using the first one.\n\n"
                fluxCalId = fluxCalId[0]

                casaCmd = casaCmd + "os.system('rm -rf %s.ampli_inf') \n"%(msName1)
                casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"
                casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                casaCmd = casaCmd + "  solint = 'inf',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                if isFullP:
                    casaCmd = casaCmd + "  refantmode = 'strict',\n"
                casaCmd = casaCmd + "  gaintype = '"+gaintypeForAmp+"',\n"
                casaCmd = casaCmd + "  calmode = 'a',\n"
                casaCmd = casaCmd + "  gaintable = ['"+bandpass+"', '"+msName1+".phase_int'])\n\n"

                if doplot == True: 
                    casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".ampli_inf', msName='"+msName1+"', interactive=False) \n\n"

                casaCmd = casaCmd + "os.system('rm -rf %s.flux_inf') \n"%(msName1)
                casaCmd = casaCmd + "os.system('rm -rf %s.fluxscale') \n"%(msName1)
                casaCmd = casaCmd + "mylogfile = casalog.logfile()\n"
                casaCmd = casaCmd + "casalog.setlogfile('"+msName1+".fluxscale')\n\n"
                casaCmd = casaCmd + fluxscaleDictName[0] + " = fluxscale(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"
                casaCmd = casaCmd + "  fluxtable = '"+msName1+".flux_inf',\n"
                casaCmd = casaCmd + "  reference = '"+str(fluxCalId)+"') # "+fieldNames[fluxCalId]+"\n\n"
                casaCmd = casaCmd + "casalog.setlogfile(mylogfile)\n\n"
                casaCmd = casaCmd + "if applyonly != True: es.fluxscale2(caltable = '"+msName1+".ampli_inf', removeOutliers=True, msName='"+msName+"', writeToFile=True, preavg=10000)\n\n"

            else:

                casaCmd = casaCmd + "os.system('rm -rf %s.flux_inf') \n"%(msName1)
                casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".flux_inf',\n"
                casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                casaCmd = casaCmd + "  solint = 'inf',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                if isFullP:
                    casaCmd = casaCmd + "  refantmode = 'strict',\n"
                casaCmd = casaCmd + "  gaintype = '"+gaintypeForAmp+"',\n"
                casaCmd = casaCmd + "  calmode = 'a',\n"
                casaCmd = casaCmd + "  gaintable = ['"+bandpass+"', '"+msName1+".phase_int'])\n\n"

                if doplot == True: casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".flux_inf', msName='"+msName1+"', interactive=False) \n\n"

        else:

            if len(fluxCalId) > 1: 
                casaCmd = casaCmd + "# Note: there is more than one Solar system object in this dataset, I'm using the first one as flux calibrator.\n\n"
            fluxCalId = fluxCalId[0]
            mytb.open(msName+'/ANTENNA')
            antList = mytb.getcol('NAME')
            mytb.close()
            print("Running es.getAntennasForFluxscale2('%s', fluxCalId='%s', refant='%s')" % (msName,str(fluxCalId),refant))
            antList1 = sfsdr.getAntennasForFluxscale2(msName, fluxCalId=str(fluxCalId), refant=refant)
            if antList1 == []:
                casalog.post("ERROR: es.getAntennasForFluxscale2 returned empty list of antennas.", 'SEVERE')
                return False


            if len(antList) == len(antList1):

                casaCmd = casaCmd + "os.system('rm -rf %s.phase_int') \n"%(msName1)
                casaCmd = casaCmd + "\ngaincal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".phase_int',\n"
                casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                casaCmd = casaCmd + "  solint = 'int',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                if isFullP:
                    casaCmd = casaCmd + "  refantmode = 'strict',\n"
                casaCmd = casaCmd + "  gaintype = 'G',\n"
                casaCmd = casaCmd + "  calmode = 'p',\n"
                casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n\n"

                if doplot == True: 
                    casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".phase_int', msName='"+msName1+"', interactive=False) \n\n"

                casaCmd = casaCmd + "os.system('rm -rf %s.ampli_inf') \n"%(msName1)
                casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"
                casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                casaCmd = casaCmd + "  solint = 'inf',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                if isFullP:
                    casaCmd = casaCmd + "  refantmode = 'strict',\n"
                casaCmd = casaCmd + "  gaintype = '"+gaintypeForAmp+"',\n"
                casaCmd = casaCmd + "  calmode = 'a',\n"
                casaCmd = casaCmd + "  gaintable = ['"+bandpass+"', '"+msName1+".phase_int'])\n\n"

                if doplot == True: 
                    casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".ampli_inf', msName='"+msName1+"', interactive=False) \n\n"

                casaCmd = casaCmd + "os.system('rm -rf %s.flux_inf') \n"%(msName1)
                casaCmd = casaCmd + "os.system('rm -rf %s.fluxscale') \n"%(msName1)
                casaCmd = casaCmd + "mylogfile = casalog.logfile()\n"
                casaCmd = casaCmd + "casalog.setlogfile('"+msName1+".fluxscale')\n\n"
                casaCmd = casaCmd + fluxscaleDictName[0] + " = fluxscale(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"
                casaCmd = casaCmd + "  fluxtable = '"+msName1+".flux_inf',\n"
                casaCmd = casaCmd + "  reference = '"+str(fluxCalId)+"') # "+fieldNames[fluxCalId]+"\n\n"
                casaCmd = casaCmd + "casalog.setlogfile(mylogfile)\n\n"
                casaCmd = casaCmd + "if applyonly != True: es.fluxscale2(caltable = '"+msName1+".ampli_inf', removeOutliers=True, msName='"+msName+"', writeToFile=True, preavg=10000)\n\n"

            else:

                if sciFieldIds[0] != '':
                    phaseCal = sfsdr.getPhaseCal(msName, valueMaps=valueMaps)  # Added by CLB
                    phaseCalNames = []  # Added by CLB
                    phaseCalIds = []
                    for i in phaseCal:  # Added by CLB
                        phaseCalNames.append(phaseCal[i]['phaseCalName']) # Added by CLB
                        phaseCalIds.append(phaseCal[i]['phaseCalId'])

                if len(antList1) < 2:
                    print('WARNING: THE SOLAR SYSTEM OBJECT SEEMS TO BE EXTREMELY RESOLVED')
                    print('WARNING: I COULD NOT FIND A SUBSET OF ANTENNAS ON WHICH TO RUN GAINCAL')
                    print('WARNING: YOU SHOULD LOOK AT THE DATA, AND THEN UPDATE THE SCRIPT')

                casaCmd += "# Note: the Solar system object used for flux calibration is highly resolved on some baselines.\n"
                casaCmd += "# Note: we will first determine the flux of the phase calibrator(s) on a subset of antennas.\n\n"

                casaCmd += "delmod"

                if sciFieldIds[0] != '':
                    casaCmd += "('%s',field='%s')\n\n" % (msName1, ",".join(map(str,list(np.unique(phaseCalIds))))) # Added by CLB
                else:
                    casaCmd += "('%s',field='%s')\n\n" % (msName1, ",".join(map(str,list(np.unique(calFieldIds)))))

                numAntList1 = len(antList1)
                antList1 = ','.join(antList1)

                casaCmd = casaCmd + "os.system('rm -rf %s.phase_short_int') \n"%(msName1)
                casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".phase_short_int',\n"
                casaCmd = casaCmd + "  field = '"+str(fluxCalId)+"', # "+fieldNames[fluxCalId]+"\n"
                casaCmd = casaCmd + "  selectdata = True,\n"
                casaCmd = casaCmd + "  antenna = '"+str(antList1)+"&',\n"
                casaCmd = casaCmd + "  solint = 'int',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                if isFullP:
                    casaCmd = casaCmd + "  refantmode = 'strict',\n"

                if numAntList1 < 5:
                    casaCmd = casaCmd + "  minblperant = "+str(numAntList1-1)+",\n"
                    casaCmd = casaCmd + "  minsnr = 2.0,\n"

                casaCmd = casaCmd + "  gaintype = 'G',\n"
                casaCmd = casaCmd + "  calmode = 'p',\n"
                casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n\n"

                ###

                calFieldIds1 = [str(i) for i in calFieldIds if i != fluxCalId]
                calFieldIds1 = ','.join(calFieldIds1)

                calFieldNames = [fieldNames[i] for i in calFieldIds if i != fluxCalId]
                calFieldNames = ','.join(calFieldNames)

                casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".phase_short_int',\n"
                casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                casaCmd = casaCmd + "  selectdata = True,\n"
                casaCmd = casaCmd + "  solint = 'int',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                if isFullP:
                    casaCmd = casaCmd + "  refantmode = 'strict',\n"

                if numAntList1 < 5:
                    casaCmd = casaCmd + "  minblperant = "+str(numAntList1-1)+",\n"
                    casaCmd = casaCmd + "  minsnr = 2.0,\n"

                casaCmd = casaCmd + "  gaintype = 'G',\n"
                casaCmd = casaCmd + "  calmode = 'p',\n"
                casaCmd = casaCmd + "  append = True,\n"
                casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n\n"

                ###

                if doplot == True: 
                    casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".phase_short_int', msName='"+msName1+"', interactive=False) \n\n"

                calFieldIds1 = [str(i) for i in calFieldIds]
                calFieldIds1 = ','.join(calFieldIds1)

                calFieldNames = [fieldNames[i] for i in calFieldIds]
                calFieldNames = ','.join(calFieldNames)

                casaCmd = casaCmd + "os.system('rm -rf %s.ampli_short_inf') \n"%(msName1)
                casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_short_inf',\n"
                casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                casaCmd = casaCmd + "  selectdata = True,\n"
                casaCmd = casaCmd + "  solint = 'inf',\n"
                casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                if isFullP:
                    casaCmd = casaCmd + "  refantmode = 'strict',\n"

                if numAntList1 < 5:
                    casaCmd = casaCmd + "  minblperant = "+str(numAntList1-1)+",\n"
                    casaCmd = casaCmd + "  minsnr = 2.0,\n"

                casaCmd = casaCmd + "  gaintype = '"+gaintypeForAmp+"',\n"
                casaCmd = casaCmd + "  calmode = 'a',\n"
                casaCmd = casaCmd + "  gaintable = ['"+bandpass+"', '"+msName1+".phase_short_int'])\n\n"

                if doplot == True: casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".ampli_short_inf', msName='"+msName1+"', interactive=False) \n\n"

                casaCmd = casaCmd + "os.system('rm -rf %s.flux_short_inf') \n"%(msName1)
                casaCmd = casaCmd + "os.system('rm -rf %s.fluxscale') \n"%(msName1)
                casaCmd = casaCmd + "mylogfile = casalog.logfile()\n"
                casaCmd = casaCmd + "casalog.setlogfile('"+msName1+".fluxscale')\n\n"
                casaCmd = casaCmd + fluxscaleDictName[0] + " = fluxscale(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_short_inf',\n"
                casaCmd = casaCmd + "  fluxtable = '"+msName1+".flux_short_inf',\n"
                casaCmd = casaCmd + "  reference = '"+str(fluxCalId)+"') # "+fieldNames[fluxCalId]+"\n\n"
                casaCmd = casaCmd + "casalog.setlogfile(mylogfile)\n\n"
                casaCmd = casaCmd + "if applyonly != True: es.fluxscale2(caltable = '"+msName1+".ampli_short_inf', removeOutliers=True, msName='"+msName+"', writeToFile=True, preavg=10000)\n\n"

                if sciFieldIds[0] != '':

                    casaCmd = casaCmd + "f = open('"+msName1+".fluxscale')\n"
                    casaCmd = casaCmd + "fc = f.readlines()\n"
                    casaCmd = casaCmd + "f.close()\n\n"

                    phaseCal = sfsdr.getPhaseCal(msName, valueMaps=valueMaps)
                    phaseCalNames = []
                    for i in phaseCal:
                        phaseCalNames.append(phaseCal[i]['phaseCalName'])

                    casaCmd = casaCmd + "for phaseCalName in "+str(list(set(phaseCalNames)))+":\n"
                    casaCmd = casaCmd + "  for i in range(len(fc)):\n"
                    casaCmd = casaCmd + "    if fc[i].find('Flux density for '+phaseCalName) != -1 and re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE) is not None:\n"

                    casaCmd = casaCmd + "      line = (re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE)).group(0)\n"
                    casaCmd = casaCmd + "      spwId = (line.split('='))[1].split()[0]\n"
                    casaCmd = casaCmd + "      flux = float((line.split(':'))[1].split()[0])\n"
                    casaCmd = casaCmd + "      setjy(vis = '"+msName1+"',\n"
                    casaCmd = casaCmd + "        field = phaseCalName.replace(';','*;').split(';')[0],\n"
                    casaCmd = casaCmd + "        spw = spwId,\n"
                    casaCmd = casaCmd + "        standard = 'manual',\n"
                    casaCmd = casaCmd + "        fluxdensity = [flux,0,0,0])\n\n"


                    casaCmd = casaCmd + "os.system('rm -rf %s.phase_int') \n"%(msName1)
                    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                    casaCmd = casaCmd + "  caltable = '"+msName1+".phase_int',\n"
                    casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                    casaCmd = casaCmd + "  solint = 'int',\n"
                    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                    if isFullP:
                        casaCmd = casaCmd + "  refantmode = 'strict',\n"
                    casaCmd = casaCmd + "  gaintype = 'G',\n"
                    casaCmd = casaCmd + "  calmode = 'p',\n"
                    casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n\n"

                    if doplot == True: 
                        casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".phase_int', msName='"+msName1+"', interactive=False) \n\n"

                    casaCmd = casaCmd + "os.system('rm -rf %s.flux_inf') \n"%(msName1)
                    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
                    casaCmd = casaCmd + "  caltable = '"+msName1+".flux_inf',\n"
                    casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
                    casaCmd = casaCmd + "  solint = 'inf',\n"
                    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
                    if isFullP:
                        casaCmd = casaCmd + "  refantmode = 'strict',\n"
                    casaCmd = casaCmd + "  gaintype = '"+gaintypeForAmp+"',\n"
                    casaCmd = casaCmd + "  calmode = 'a',\n"
                    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"', '"+msName1+".phase_int'])\n\n"

                    if doplot == True: 
                        casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".flux_inf', msName='"+msName1+"', interactive=False) \n\n"


        if sciFieldIds[0] != '' and calmode2 == 'ap':
            casaCmd = casaCmd + "os.system('rm -rf %s.phase_inf') \n"%(msName1)
            casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  caltable = '"+msName1+".phase_inf',\n"
            casaCmd = casaCmd + "  field = '"+calFieldIds1+"', # "+calFieldNames+"\n"
            casaCmd = casaCmd + "  solint = 'inf',\n"
            casaCmd = casaCmd + "  refant = '"+refant+"',\n"
            if isFullP:
                casaCmd = casaCmd + "  refantmode = 'strict',\n"
            casaCmd = casaCmd + "  gaintype = 'G',\n"
            casaCmd = casaCmd + "  calmode = 'p',\n"
            casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n\n"

            if doplot == True: 
                casaCmd = casaCmd + "if applyonly != True: es.checkCalTable('"+msName1+".phase_inf', msName='"+msName1+"', interactive=False) \n"

    return casaCmd

#####################################
def doB2BGainCalibrationPartI(msName, msName1='', refant='', bandpass='', doplot=True, iHaveSplitMyScienceSpw=False, valueMaps={},
                              combineB2BLFspws=False, combineB2BLFHFspws=False, combineB2BDGCspws=False):
    """Generate code for the first gain calibration step of a calibration script for B2B data."""

    print('\n*** doB2BGainCalibrationPartI ***')
    print('Gathering information ...')

    if msName1 == '': msName1 = msName
    if refant == '': 
        casalog.post('ERROR: No reference antenna specified.', 'SEVERE')
        return False
    if bandpass == '': 
        casalog.post('ERROR: No bandpass cal table specified.', 'SEVERE')
        return False

    ### determine the field, scan, and spw ids to be used in the code

    if msName in valueMaps.keys():
        vm = valueMaps[msName]
        print('Using canned ValueMap.')
    else:
        vm = aU.ValueMapping(msName)
        valueMaps[msName] = vm

    intentSources = sfsdr.getIntentsAndSourceNames(msName)

    diffGainCalId = intentSources['CALIBRATE_DIFFGAIN']['id'][0]

    spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#SIGNAL')
    if spwsDiffgainsig == []:
        spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#ON_SOURCE')
        sigintent = 'CALIBRATE_DIFFGAIN#ON_SOURCE'
    else:
        sigintent = 'CALIBRATE_DIFFGAIN#SIGNAL'

    mymsmd = msmdtool()
    mymsmd.open(msName)

    diffGainCalScanListLow  = ','.join([str(i) for i in mymsmd.scansforintent('CALIBRATE_DIFFGAIN#REFERENCE')])
    diffGainCalScanListHigh = ','.join([str(i) for i in mymsmd.scansforintent(sigintent)])

    # added all diffgainscals - ATM are alone so never have DIFFGAIN intent 
    diffGainCalScanList = ','.join([str(i) for i in mymsmd.scansforintent('CALIBRATE_DIFFGAIN#*')])
    
    mymsmd.close()

    spwHighDict = sfsdr.getSpwInfo(msName,intent=sigintent, caching=True)
    spwLowDict = sfsdr.getSpwInfo(msName,intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
    spwHigh = sorted(spwHighDict.keys())
    spwLow = sorted(spwLowDict.keys())

    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/OBSERVATION')
    obsTimeRange = mytb.getcol('TIME_RANGE')
    mytb.close()
    obsTimeStart = ((obsTimeRange[0]/86400.0)+2400000.5-2440587.5)*86400.0
    obsTimeStart = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTimeStart[0])))

    # determine spw mapping
    # conforming to correct regimes with extra mapping and improved naming
    
    if iHaveSplitMyScienceSpw == True:
        spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
        spwIds = sorted(spwInfo.keys())

        spwHighSplit = []
        for myspw in spwHigh:
            spwHighSplit.append(spwIds.index(myspw))
        spwHighStr = ','.join([str(i) for i in spwHighSplit])

        spwLowSplit = []
        for myspw in spwLow:
            spwLowSplit.append(spwIds.index(myspw))
        spwLowStr = ','.join([str(i) for i in spwLowSplit])


        # need to now consider the cases for mapping

        # logic is simplified as we don't need all the maps defined here
        # we need ONLY consider LFtoHF B2B transfer
        # and HFtoHF_DGC in Case D
        # trigger is the same on combineB2BLFspws = True
        # and now combineB2BDGCspws = True
        
        LFtoHF = list(range(max(spwHighSplit)+1))   # only need this for B2B transfer

            
        lowestSpwLow = int(min(spwLowSplit))
        
        if obsTimeStart > '2021-03-01T00:00:00':
            for spwRep in spwHighSplit: 
                BBnum = spwHighDict[spwRep]['basebandNum'] 
                for spwMap in spwLowSplit:  
                    if spwLowDict[spwMap]['basebandNum']==BBnum:
                        if combineB2BLFspws:
                            LFtoHF[spwRep] = lowestSpwLow
                        else:
                            LFtoHF[spwRep] = spwMap  

        else:
            for spwId,spwRep in enumerate(spwHighSplit):
                if combineB2BLFspws:
                    LFtoHF[spwRep] = lowestSpwLow
                else:
                    LFtoHF[spwRep] = spwLowSplit[spwId]

        # Now Case D
        # these are HF to HF self
        if combineB2BDGCspws:
            HFtoHF_DGC = list(range(max(spwHighSplit)+1))
            for spwMap in spwHighSplit:
                HFtoHF_DGC[spwMap] = min(spwHighSplit)  # maps to min of HF index



                    
    else: # no reindexing took place
        spwHighStr = ','.join([str(i) for i in spwHigh])
        spwLowStr = ','.join([str(i) for i in spwLow])

        LFtoHF = list(range(max(spwHigh)+1))
        lowestSpwLow = int(min(spwLow))

        if obsTimeStart > '2021-03-01T00:00:00':
            for spwRep in spwHigh: 
                BBnum = spwHighDict[spwRep]['basebandNum']
                for spwMap in spwLow:  
                    if spwLowDict[spwMap]['basebandNum']==BBnum:
                        if combineB2BLFspws:
                            LFtoHF[spwRep] = lowestSpwLow
                        else:
                            LFtoHF[spwRep] = spwMap
        else:
            if combineB2BLFspws: 
                for spwRep in spwHigh:
                    LFtoHF[spwRep] = lowestSpwLow
            else:
                for spwId,spwRep in enumerate(spwHigh):
                    LFtoHF[spwRep] = spwLow[spwId]

        # Now Case D
        # these are HF to HF self
        if combineB2BDGCspws:
            HFtoHF_DGC = list(range(max(spwHigh)+1))
            for spwMap in spwHigh:
                HFtoHF_DGC[spwMap] = min(spwHigh)  # maps to min of HF index


                    
    ###

    print('diffGainCalScanListLow ', diffGainCalScanListLow)
    print('diffGainCalScanListHigh ', diffGainCalScanListHigh)
    print('spwHighStr ', spwHighStr)
    print('spwLowStr ', spwLowStr)
    print('LFtoHF ', str(LFtoHF))

    if combineB2BDGCspws:
        print('Combine HF SpWs for B2B solution')
        print('HFtoHF_DGC ', str(HFtoHF_DGC))

    ###

    print('Writing code ...')

    casaCmd = ''

    casaCmd = casaCmd + "# do a phase offset so all respective HF and LF SPW can be merged for better SNR, if required\n" 
    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".DGCphase_phasediff')\n"
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".DGCphase_phasediff',\n"
    casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"   
    casaCmd = casaCmd + "  scan = '"+diffGainCalScanList+"',\n"
    casaCmd = casaCmd + "  solint = 'inf',\n"
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  combine = 'scan',\n"  
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    casaCmd = casaCmd + "  gaintable = '"+bandpass+"')\n"

    if doplot:
        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".DGCphase_phasediff', msName='"+msName1+"', interactive=False)\n" 


    casaCmd = casaCmd + "\n# LF temporal solutions fast\n"
    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".DGCphaselow_int')\n" 

    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".DGCphaselow_int',\n"
    casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
    casaCmd = casaCmd + "  scan = '"+diffGainCalScanListLow+"',\n"
    casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
    casaCmd = casaCmd + "  solint = 'int',\n"
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    if combineB2BLFspws:
        casaCmd = casaCmd + "  combine = 'spw',\n"
    casaCmd = casaCmd + "  refantmode = 'strict',\n"
    casaCmd = casaCmd + "  interp = ['linear','linear'],\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"', '"+msName1+".DGCphase_phasediff'],\n"
    casaCmd = casaCmd + "  gainfield = ['',''])\n"

    if doplot:
        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".DGCphaselow_int', msName='"+msName1+"', interactive=False)\n" 

    casaCmd = casaCmd + "\n# LF temporal solutions slow\n"
    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".DGCphaselow_inf')\n"
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".DGCphaselow_inf',\n"
    casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
    casaCmd = casaCmd + "  scan = '"+diffGainCalScanListLow+"',\n"
    casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
    casaCmd = casaCmd + "  solint = 'inf',\n"
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    if combineB2BLFspws:
        casaCmd = casaCmd + "  combine = 'spw',\n"
    casaCmd = casaCmd + "  refantmode = 'strict',\n"
    casaCmd = casaCmd + "  interp = ['linear','linear'],\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff'],\n"
    casaCmd = casaCmd + "  gainfield = ['',''])\n"

    if doplot:
        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".DGCphaselow_inf', msName='"+msName1+"', interactive=False)\n"


    casaCmd = casaCmd + "\n## HF fast solutions\n"

    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".DGCphasehigh_int')\n"

    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".DGCphasehigh_int',\n"
    casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
    casaCmd = casaCmd + "  scan= '"+diffGainCalScanListHigh+"',\n"  
    casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
    casaCmd = casaCmd + "  solint = 'int',\n"
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    if combineB2BLFHFspws:
        casaCmd = casaCmd + "  combine = 'spw',\n"
    casaCmd = casaCmd + "  refantmode = 'strict',\n"
    casaCmd = casaCmd + "  interp = ['linear','linear'],\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff'],\n"
    casaCmd = casaCmd + "  gainfield = ['',''])\n"


    if doplot:
        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".DGCphasehigh_int', msName='"+msName1+"', interactive=False)\n" 


    casaCmd = casaCmd + "\n# get the DGC band-offset\n"
    casaCmd = casaCmd + "# apply the LF to the HF and find the offset for each DGC 'group'\n"
    casaCmd = casaCmd + "# if there are multiple DGC sequences - this does each one independently (hence 'multi')\n"

    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".DGCoffset_multi')\n"
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".DGCoffset_multi',\n"
    casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
    casaCmd = casaCmd + "  scan = '"+diffGainCalScanListHigh+"',\n"  
    casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
    casaCmd = casaCmd + "  solint = '8min',  # max. duration of one DGC visit block\n"
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  refantmode = 'strict',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  minsnr=3.0,\n"
    if combineB2BDGCspws:
        casaCmd = casaCmd + "  combine = 'scan,spw',\n"
    else:
        casaCmd = casaCmd + "  combine='scan',\n"
    casaCmd = casaCmd + "  interp=['linear','linear','linearPD'],\n"
    casaCmd = casaCmd + "  spwmap=[[],[],"+str(LFtoHF)+"],\n"
    casaCmd = casaCmd + "  gainfield = ['','',''],\n"  
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".DGCphaselow_inf'])\n"

    if doplot:
        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".DGCoffset_multi', msName='"+msName1+"', interactive=False)\n"

        casaCmd = casaCmd + "\nif applyonly != True:\n"
        casaCmd = casaCmd + "  # QA2 check to verify the slow (inf) LF to HF and offset worked\n"
        casaCmd = casaCmd + "  # apply the LF to the HF and the multi offset\n"
        casaCmd = casaCmd + "  os.system('rm -rf "+msName1+".DGCresidual_offset_multi')\n"
        casaCmd = casaCmd + "  gaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "    caltable = '"+msName1+".DGCresidual_offset_multi',\n"
        casaCmd = casaCmd + "    field = '"+str(diffGainCalId)+"',\n"
        casaCmd = casaCmd + "    scan = '"+diffGainCalScanListHigh+"',\n"
        casaCmd = casaCmd + "    spw = '"+spwHighStr+"',\n"
        casaCmd = casaCmd + "    solint = 'inf',\n"
        casaCmd = casaCmd + "    refant = '"+refant+"',\n"
        casaCmd = casaCmd + "    refantmode = 'strict',\n"
        casaCmd = casaCmd + "    gaintype = 'G',\n"
        casaCmd = casaCmd + "    calmode = 'p',\n"
        if combineB2BLFHFspws:
            casaCmd = casaCmd + "    #combine='spw', # NOTE: uncomment if SNR is too low for the residual solns to be useful in checking\n"
        casaCmd = casaCmd + "    minsnr=3.0,\n"
        casaCmd = casaCmd + "    interp=['linear','linear','linearPD','linear'],\n"
        if combineB2BDGCspws:
            casaCmd = casaCmd + "    spwmap=[[],[],"+str(LFtoHF)+","+str(HFtoHF_DGC)+"],\n"
        else:
            casaCmd = casaCmd + "    spwmap=[[],[],"+str(LFtoHF)+",[]],\n"
        casaCmd = casaCmd + "    gainfield = ['','','',''],\n"
        casaCmd = casaCmd + "    gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".DGCphaselow_inf','"+msName1+".DGCoffset_multi'])\n"

        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".DGCresidual_offset_multi', msName='"+msName1+"', interactive=False)\n"

        casaCmd = casaCmd + "\nif applyonly != True:\n"
        casaCmd = casaCmd + "  # QA2 check to verify the fast (int) LF to HF and offset worked\n"
        casaCmd = casaCmd + "  os.system('rm -rf "+msName1+".DGCresidual_offset_multi_int')\n" 
        casaCmd = casaCmd + "  gaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "    caltable = '"+msName1+".DGCresidual_offset_multi_int',\n"
        casaCmd = casaCmd + "    field = '"+str(diffGainCalId)+"',\n"
        casaCmd = casaCmd + "    scan = '"+diffGainCalScanListHigh+"',\n"
        casaCmd = casaCmd + "    spw = '"+spwHighStr+"',\n"
        casaCmd = casaCmd + "    solint = 'int',\n"
        casaCmd = casaCmd + "    refant = '"+refant+"',\n"
        casaCmd = casaCmd + "    refantmode = 'strict',\n"
        casaCmd = casaCmd + "    gaintype = 'G',\n"
        casaCmd = casaCmd + "    calmode = 'p',\n"
        if combineB2BLFHFspws:
            casaCmd = casaCmd + "    combine='spw', # SNR is probably too low; need combine for 'int' checking\n"        
        casaCmd = casaCmd + "    minsnr=3.0,\n"
        casaCmd = casaCmd + "    interp=['linear','linear','linearPD','linear'],\n"
        if combineB2BDGCspws:
            casaCmd = casaCmd + "    spwmap=[[],[],"+str(LFtoHF)+","+str(HFtoHF_DGC)+"],\n"
        else:
            casaCmd = casaCmd + "    spwmap=[[],[],"+str(LFtoHF)+",[]],\n"
        casaCmd = casaCmd + "    gainfield = ['','','',''],\n"
        casaCmd = casaCmd + "    gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".DGCphaselow_inf','"+msName1+".DGCoffset_multi'])\n"

        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".DGCresidual_offset_multi_int', msName='"+msName1+"', interactive=False)\n"
    #endif doplot


    return casaCmd

#####################################
def doB2BGainCalibrationPartII(msName, msName1='', refant='', bandpass='', doplot=True, ampForSci=[], iHaveSplitMyScienceSpw=False, valueMaps={},
                               combineB2BLFspws=False, combineB2BLFHFspws=False, combineB2BDGCspws=False):
    """Generate code for the second gain calibration step of a calibration script for B2B data."""

    print('\n*** doB2BGainCalibrationPartII ***')
    print('Gathering information ...')

    if msName1 == '': msName1 = msName
    if refant == '': 
        casalog.post('ERROR: No reference antenna specified.', 'SEVERE')
        return False
    if bandpass == '': 
        casalog.post('ERROR: No bandpass cal table specified.', 'SEVERE')
        return False

    ### determine the scan lists and the field and spw ids to be used in the code

    if msName in valueMaps.keys():
        vm = valueMaps[msName]
        print('Using canned ValueMap.')
    else:
        vm = aU.ValueMapping(msName)
        valueMaps[msName] = vm

    intentSources = sfsdr.getIntentsAndSourceNames(msName)

    diffGainCalId = intentSources['CALIBRATE_DIFFGAIN']['id'][0]
    phaseCalId  = intentSources['CALIBRATE_PHASE']['id'][0]
    i = 0
    while phaseCalId == diffGainCalId and i < len(intentSources['CALIBRATE_PHASE']['id'])-1:
        i += 1
        phaseCalId  = intentSources['CALIBRATE_PHASE']['id'][i]

    fluxCalId  = intentSources['CALIBRATE_FLUX']['id'][0]

    fcCoversHFandLF = False
    if intentSources['CALIBRATE_FLUX']['spw'] == intentSources['CALIBRATE_DIFFGAIN']['spw']: # fluxcal covers HF and LF
        fcCoversHFandLF = True
        fluxCalLFId = fluxCalId
    else:
        fluxCalLFId = diffGainCalId 

    bpCalId  = intentSources['CALIBRATE_BANDPASS']['id'][0]

    fieldCal = set((bpCalId, diffGainCalId, phaseCalId))
    fieldCalStr = ','.join([str(i) for i in fieldCal])
    fieldCal2 = set((fluxCalId, bpCalId, diffGainCalId, phaseCalId))
    fieldCal2Str = ','.join([str(i) for i in fieldCal2])
    fieldCal3 = set((fluxCalId, bpCalId, diffGainCalId))
    fieldCal3Str = ','.join([str(i) for i in fieldCal3])


    spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#SIGNAL')
    if spwsDiffgainsig == []:
        spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#ON_SOURCE')
        sigintent = 'CALIBRATE_DIFFGAIN#ON_SOURCE'
    else:
        sigintent = 'CALIBRATE_DIFFGAIN#SIGNAL'

    mymsmd = msmdtool()
    mymsmd.open(msName)
    fluxCalScanList = ','.join([str(i) for i in mymsmd.scansforintent('CALIBRATE_FLUX#ON_SOURCE')])
    diffGainCalScansLow = [i for i in mymsmd.scansforintent('CALIBRATE_DIFFGAIN#REFERENCE')]
    diffGainCalScansHigh = [i for i in mymsmd.scansforintent(sigintent)]
    diffGainCalScanListLow_groups = [','.join(map(str, group))
                                     for group in np.split(diffGainCalScansLow,
                                                           np.where(np.diff(diffGainCalScansLow) > 3)[0] + 1)]
    diffGainCalScanListHigh_groups = [','.join(map(str, group))
                                     for group in np.split(diffGainCalScansHigh,
                                                           np.where(np.diff(diffGainCalScansHigh) > 3)[0] + 1)]
    bandpassScans = [str(i) for i in mymsmd.scansforintent('CALIBRATE_BANDPASS#ON_SOURCE')]
    mymsmd.close()
    if len(bandpassScans) != 2:
        casalog.post('ERROR: Found '+str(len(bandpassScans))+' BANDPASS scans:'+str(bandpassScans)+', expected two (one LF, one HF).', 'SEVERE')
        return False
    else:
        bandpassScanListLow = bandpassScans[0]
        bandpassScanListHigh = bandpassScans[1]
    
    if fluxCalLFId == fluxCalId and fcCoversHFandLF:
        fluxCalLFScanList = fluxCalScanList
    else:
        fluxCalLFScanList = diffGainCalScanListLow

    spwHighDict = sfsdr.getSpwInfo(msName,intent=sigintent, caching=True)
    spwLowDict = sfsdr.getSpwInfo(msName,intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
    spwHigh = sorted(spwHighDict.keys())
    spwLow = sorted(spwLowDict.keys())

    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/OBSERVATION')
    obsTimeRange = mytb.getcol('TIME_RANGE')
    mytb.close()
    obsTimeStart = ((obsTimeRange[0]/86400.0)+2400000.5-2440587.5)*86400.0
    obsTimeStart = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTimeStart[0])))

    ## only doing LF combine 
    ## for accounting for possible HF combine that can be done via CASE C
    
    if iHaveSplitMyScienceSpw == True:
        spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
        spwIds = sorted(spwInfo.keys())

        spwHighSplit = []
        for myspw in spwHigh:
            spwHighSplit.append(spwIds.index(myspw))
        spwHighStr = ','.join([str(i) for i in spwHighSplit])

        spwLowSplit = []
        for myspw in spwLow:
            spwLowSplit.append(spwIds.index(myspw))
        spwLowStr = ','.join([str(i) for i in spwLowSplit])
        
        LFtoLF = list(range(max(spwLowSplit)+1))
        HFtoHF = list(range(max(spwHighSplit)+1))


        # no baseband mapping required for LFtoLF and HftoHF
        if combineB2BLFspws:
            for spwMap in spwLowSplit:
                LFtoLF[spwMap] =  int(min(spwLowSplit))      
        if combineB2BLFHFspws:
            for spwMap in spwHighSplit:
                HFtoHF[spwMap] =  int(min(spwHighSplit))                                       

                
    else: # no reindexing took place
        spwHighStr = ','.join([str(i) for i in spwHigh])
        spwLowStr = ','.join([str(i) for i in spwLow])

        LFtoLF = list(range(max(spwLow)+1))
        HFtoHF = list(range(max(spwHigh)+1))                                      

        # no baseband mapping required for LFtoLF and HFtoHF
        if combineB2BLFspws:
            for spwMap in spwLow:
                LFtoLF[spwMap] =  int(min(spwLow))     
        if combineB2BLFHFspws:
            for spwMap in spwHigh:
                HFtoHF[spwMap] =  int(min(spwHigh))                   


    ###


    print('fluxCalScanList ', fluxCalScanList)
    print('fluxCalLFScanList ', fluxCalLFScanList)
    print('spwHighStr ', spwHighStr)
    print('spwLowStr ', spwLowStr)

    if combineB2BLFspws:
        print('LFtoLF ', str(LFtoLF))
    if combineB2BLFHFspws:
        print('HFtoHF ', str(HFtoHF))                                      
    ###

    print('Writing code ...')

    casaCmd = ''

    casaCmd = casaCmd + "# LF fast phase as a QA2 check for stability scan drop-outs\n"
    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".phaselow_int')\n" 
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".phaselow_int',\n" 
    casaCmd = casaCmd + "  field = '"+fieldCalStr+"',\n"
    casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
    casaCmd = casaCmd + "  solint = 'int',\n"
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  refantmode = 'strict',\n"
    if combineB2BLFspws:
        casaCmd = casaCmd + "  combine='spw',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff'])\n"

    casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".phaselow_int', msName='"+msName1+"', interactive=False)\n" 

    casaCmd = casaCmd + "\n# LF slow phase solution\n"
    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".phaselow_inf')\n" 
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".phaselow_inf',\n"
    casaCmd = casaCmd + "  field = '"+fieldCalStr+"',\n"    
    casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
    casaCmd = casaCmd + "  solint = 'inf',\n"  
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  refantmode='strict',\n"
    if combineB2BLFspws:
        casaCmd = casaCmd + "  combine='spw',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff'])\n"
 
    casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".phaselow_inf', msName='"+msName1+"', interactive=False)\n"

    casaCmd = casaCmd + "\n## HF phase up on all strong sources, i.e. the BP, flux cal, and DGC for amp gain table\n"
    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".phasehigh_int')\n" 
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".phasehigh_int',\n" 
    casaCmd = casaCmd + "  field = '"+fieldCal3Str+"',\n"    
    casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
    casaCmd = casaCmd + "  solint = 'int',\n"  
    casaCmd = casaCmd + "  refant = '"+refant+"',\n"
    casaCmd = casaCmd + "  refantmode = 'strict',\n"
    casaCmd = casaCmd + "  gaintype = 'G',\n"
    casaCmd = casaCmd + "  calmode = 'p',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    if combineB2BLFHFspws:
        casaCmd = casaCmd + "  combine='spw',\n"                                      
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff'])\n"
        
    casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".phasehigh_int', msName='"+msName1+"', interactive=False)\n" 

    casaCmd = casaCmd + "\n# Amplitude calibration HF\n"
    casaCmd = casaCmd + "os.system('rm -rf "+msName1+".ampli_inf')\n" 
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"  
    if fluxCalId == bpCalId:
        casaCmd = casaCmd + "  field = '"+str(fluxCalId)+"',\n"
        casaCmd = casaCmd + "  scan = '"+bandpassScanListHigh+"',\n" 
    else:
        casaCmd = casaCmd + "  field = '"+str(fluxCalId)+","+str(bpCalId)+"',\n"
        casaCmd = casaCmd + "  scan = '"+fluxCalScanList+","+bandpassScanListHigh+"',\n" 
    casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
    casaCmd = casaCmd + "  solint = 'inf',\n"  
    casaCmd = casaCmd + "  combine = 'scan',\n"  
    casaCmd = casaCmd + "  refant = '"+refant+"',\n" 
    casaCmd = casaCmd + "  gaintype = 'T',\n"
    casaCmd = casaCmd + "  calmode = 'a',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phasehigh_int'],\n"
    if combineB2BLFHFspws:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(HFtoHF)+"],\n"
    casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','',''])\n"

    for i in range(len(diffGainCalScanListHigh_groups)):
        casaCmd = casaCmd + "\n# Amplitude calibration HF - diffgain, append\n"
        casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"  
        casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
        casaCmd = casaCmd + "  scan = '"+diffGainCalScanListHigh_groups[i]+"',\n" 
        casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
        casaCmd = casaCmd + "  solint = 'inf',\n"  
        casaCmd = casaCmd + "  combine = 'scan',\n"  
        casaCmd = casaCmd + "  refant = '"+refant+"',\n" 
        casaCmd = casaCmd + "  gaintype = 'T',\n"
        casaCmd = casaCmd + "  calmode = 'a',\n"
        casaCmd = casaCmd + "  minsnr = 3.0,\n"
        casaCmd = casaCmd + "  append = True,\n"
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phasehigh_int'],\n"
        if combineB2BLFHFspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(HFtoHF)+"],\n"
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','',''])\n"

    
    casaCmd = casaCmd + "\n# Amplitude calibration LF, append\n"
    casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"  
    if fluxCalId == bpCalId:
        casaCmd = casaCmd + "  field = '"+str(fluxCalId)+"',\n"
        casaCmd = casaCmd + "  scan = '"+bandpassScanListLow+"',\n" 
    else:
        casaCmd = casaCmd + "  field = '"+str(fluxCalId)+","+str(bpCalId)+"',\n"
        casaCmd = casaCmd + "  scan = '"+fluxCalScanList+","+bandpassScanListLow+"',\n" 
    casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
    casaCmd = casaCmd + "  solint = 'inf',\n"  
    casaCmd = casaCmd + "  combine = 'scan',\n"  
    casaCmd = casaCmd + "  refant = '"+refant+"',\n" 
    casaCmd = casaCmd + "  gaintype = 'T',\n"
    casaCmd = casaCmd + "  calmode = 'a',\n"
    casaCmd = casaCmd + "  minsnr = 3.0,\n"
    casaCmd = casaCmd + "  append = True,\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_int'],\n"
    if combineB2BLFspws:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoLF)+"],\n"
    casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','',''])\n"

    for i in range(len(diffGainCalScanListLow_groups)):
        casaCmd = casaCmd + "\n# Amplitude calibration LF - diffgain, append\n"
        casaCmd = casaCmd + "gaincal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"  
        casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
        casaCmd = casaCmd + "  scan = '"+diffGainCalScanListLow_groups[i]+"',\n" 
        casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
        casaCmd = casaCmd + "  solint = 'inf',\n"  
        casaCmd = casaCmd + "  combine = 'scan',\n"  
        casaCmd = casaCmd + "  refant = '"+refant+"',\n" 
        casaCmd = casaCmd + "  gaintype = 'T',\n"
        casaCmd = casaCmd + "  calmode = 'a',\n"
        casaCmd = casaCmd + "  minsnr = 3.0,\n"
        casaCmd = casaCmd + "  append = True,\n"
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_int'],\n"
        if combineB2BLFspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoLF)+"],\n"
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','',''])\n"

    
    casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+msName1+".ampli_inf', msName='"+msName1+"', interactive=False)\n" 
    if fluxCalId != diffGainCalId:
        casaCmd = casaCmd + "\nos.system('rm -rf %s.flux_inf') \n"%(msName1)
        casaCmd = casaCmd + "os.system('rm -rf %s.fluxscale') \n"%(msName1)
        casaCmd = casaCmd + "mylogfile = casalog.logfile()\n"
        casaCmd = casaCmd + "casalog.setlogfile('"+msName1+".fluxscale')\n\n"
        casaCmd = casaCmd + "fluxscaleDict = fluxscale(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msName1+".ampli_inf',\n"
        casaCmd = casaCmd + "  fluxtable = '"+msName1+".flux_inf',\n"
        casaCmd = casaCmd + "  reference = '"+str(fluxCalId)+"') # fluxcal\n\n"
        casaCmd = casaCmd + "casalog.setlogfile(mylogfile)\n\n"


    return casaCmd

#####################################

def doApplyBandpassAndGainCalTables(msName, msName1='', iHaveSplitMyScienceSpw=False, bandpass='', phaseForCal='', phaseForSci='', flux='', phaseDiffCalTableName='', ampForSci='', useForLoop=True, renorm='', valueMaps={}):
    """Generate code for the final applycal step for the given MS (standard cal setup and BWSW case)"""

    print('\n*** doApplyBandpassAndGainCalTables ***')
    print('Gathering information ...')

    casaCmd = ''

    if bandpass == '' or phaseForCal == '' or phaseForSci == '' or flux == '': sys.exit('ERROR: Missing table(s).')
    if msName1 == '': msName1 = msName
    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/FIELD')
    fieldNames = mytb.getcol('NAME')
    mytb.close()

    fieldIds = list(range(len(fieldNames)))

    intentSources = sfsdr.getIntentsAndSourceNames(msName)
    sciFieldIds = intentSources['OBSERVE_TARGET']['id']
    if sciFieldIds[0] == '': sys.exit('ERROR: There seems to be no science field.')

    checkFieldIds = intentSources['OBSERVE_CHECK']['id']
    if checkFieldIds[0] != '':
        sciFieldIds += intentSources['OBSERVE_CHECK']['id']

    calFieldIds = [i for i in fieldIds if i not in sciFieldIds]
    if len(calFieldIds) == '': sys.exit('ERROR: There seems to be no calibrator field.')

    mymsmd = msmdtool()
    mymsmd.open(msName)
    hasdata = []
    for i in calFieldIds:
        calFieldIntents = mymsmd.intentsforfield(i)
        hasdata1 = 0
        for j in calFieldIntents:
            if re.search('^CALIBRATE_(POINTING|ATMOSPHERE|WVR)', j) == None:
                hasdata1 = 1
                break
        hasdata.append(hasdata1)
    mymsmd.close()

    calFieldIds = [calFieldIds[i] for i in range(len(calFieldIds)) if hasdata[i] == 1]

    if len(phaseDiffCalTableName) != 0:

        # calspwmap creation for spw-combined and non-spw-combined phase_int

        mymsmd = msmdtool()
        mymsmd.open(msName)

        spwInfo3 = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
        spwIds3 = sorted(spwInfo3.keys())
        spwInfoDGCR = sfsdr.getSpwInfo(msName, intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
        spwIdsDGCR = sorted(spwInfoDGCR.keys())
        spwInfoScience = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET#ON_SOURCE', caching=True)
        spwIdsScience = sorted(spwInfoScience.keys())

        if iHaveSplitMyScienceSpw == True:
            numSpws = len(spwIds3)
            tmpIds = []
            for i in spwIdsDGCR:
                tmpIds.append(spwIds3.index(i))
            spwIdsDGCR = sorted(tmpIds)
            tmpIds = []
            for i in spwIdsScience:
                tmpIds.append(spwIds3.index(i))
            spwIdsScience = sorted(tmpIds)
            tmpIds = []
            for i in spwIds3:
                tmpIds.append(spwIds3.index(i))
            spwIds3 = sorted(tmpIds)

        else:
            numSpws = max(spwIds3)+2

        calspwmapcomb = [] # version to be used with spw-combination

        # decide whether the narrow or the wide SPWs come first
        narrowFirst = True
        for i in spwIdsScience:
            for j in spwIdsDGCR: 
                if i > j:
                    narrowFirst = False
                    break

        # create calspwmapcomb accordingly
        minDGCRSpw = min(spwIdsDGCR)
        minSciSpw = min(spwIdsScience)
        if narrowFirst:
            for i in range(minDGCRSpw):
                calspwmapcomb.append(minSciSpw)
            for i in range(numSpws-minDGCRSpw):
                calspwmapcomb.append(minDGCRSpw)
        else:
            for i in range(minSciSpw):
                calspwmapcomb.append(minDGCRSpw)
            for i in range(numSpws-minSciSpw):
                calspwmapcomb.append(minSciSpw)


        ## create calspwmap for phase and flux 
        calspwmapf = []

        for i in range(numSpws):
            if i in spwIdsScience:
                myDGCRSpw = spwIdsDGCR[mymsmd.baseband(i)-1] # wide SPWs are ordered by baseband number
                calspwmapf.append(myDGCRSpw)
            else:
                calspwmapf.append(i)

        mymsmd.close()        


    phaseCal = sfsdr.getPhaseCal(msName, valueMaps=valueMaps)
    phaseCalFieldIds = []
    for j in phaseCal:
        phaseCalFieldIds.append(phaseCal[j]['phaseCalId'])
    phaseCalFieldIds = sorted(dict.fromkeys(phaseCalFieldIds).keys())
    calFieldIds = [i for i in calFieldIds if i not in phaseCalFieldIds]
    calFieldIds1 = [str(i) for i in calFieldIds]

    calFieldNames = [fieldNames[i] for i in calFieldIds]
    calFieldNames = ','.join(calFieldNames)

    gainTable = []
    gainTable.append(bandpass)

    if len(phaseDiffCalTableName) == 0:

        gainTable.append(phaseForCal)
        gainTable.append(flux)

        if useForLoop == True:
            casaCmd = casaCmd + "for i in "+str(calFieldIds1)+": # "+calFieldNames+"\n"
            casaCmd = casaCmd + "  applycal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "    field = str(i),\n"
            casaCmd = casaCmd + "    gaintable = "+str(gainTable)+",\n"
            casaCmd = casaCmd + "    gainfield = ['', i, i],\n"
            if re.search('^3.3', aU.getCasaVersion()) == None:
                casaCmd = casaCmd + "    interp = 'linear,linear',\n"
            else:
                casaCmd = casaCmd + "    interp = 'linear',\n"

            casaCmd = casaCmd + "    calwt = True,\n"
            casaCmd = casaCmd + "    flagbackup = False)\n"
        else:
            for i in calFieldIds:
                casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  field = '"+str(i)+"', # "+fieldNames[i]+"\n"
                casaCmd = casaCmd + "  gaintable = "+str(gainTable)+",\n"
                casaCmd = casaCmd + "  gainfield = ['', '"+str(i)+"', '"+str(i)+"'],\n"
                if re.search('^3.3', aU.getCasaVersion()) == None:
                    casaCmd = casaCmd + "  interp = 'linear,linear',\n"
                else:
                    casaCmd = casaCmd + "  interp = 'linear',\n"
            
                casaCmd = casaCmd + "  calwt = True,\n"
                casaCmd = casaCmd + "  flagbackup = False)\n\n"

    else:

            gainTable.append(phaseDiffCalTableName[0])

            bpassCalId = intentSources['CALIBRATE_BANDPASS']['id']
            phaseCalId = intentSources['CALIBRATE_PHASE']['id']
            ampCalId = intentSources['CALIBRATE_AMPLI']['id'] + intentSources['CALIBRATE_FLUX']['id']
            ampCalId = [i for i in ampCalId if i != '']

            if bpassCalId[0] != '':
                if len(bpassCalId) != 1: casaCmd = casaCmd + "# Note: there is more than one bandpass calibrator, I'm picking the first one: "+fieldNames[bpassCalId[0]]+".\n"
                bpassCalId = bpassCalId[0]

            else:
                casaCmd = casaCmd + "# Note: there is no bandpass calibrator, I'm picking a phase calibrator.\n"
                phaseOnlyCalId = [i for i in phaseCalId if i not in ampCalId]
                if len(phaseOnlyCalId) != 1: casaCmd = casaCmd + "# Note: there is more than one phase calibrator, I'm picking the first one: "+fieldNames[phaseOnlyCalId[0]]+".\n"
                bpassCalId = phaseOnlyCalId[0]


            casaCmd = casaCmd + "# Case combine='spw' in phase_int\n"
            casaCmd = casaCmd + "#calspwmap = {"+str(bpassCalId)+": "+repr(calspwmapcomb)
            for i in calFieldIds + phaseCalFieldIds:
                if i == bpassCalId: continue
                if i in phaseCalFieldIds:
                    casaCmd = casaCmd + ",\n" + "#             "+str(i)+": "+repr([minDGCRSpw for i in range(numSpws)])
                else:
                    casaCmd = casaCmd + ",\n" + "#             "+str(i)+": "+repr(calspwmapcomb)
            casaCmd = casaCmd + "}\n\n"

            casaCmd = casaCmd + "# Case combine='' in phase_int\n"
            casaCmd = casaCmd + "calspwmap = {"+str(bpassCalId)+": list(range("+str(numSpws)+"))"
            for i in calFieldIds+phaseCalFieldIds:
                if i == bpassCalId: continue
                if i in phaseCalFieldIds:
                    casaCmd = casaCmd + ",\n" + "             "+str(i)+": "+repr(calspwmapf)
                else:
                    casaCmd = casaCmd + ",\n" + "             "+str(i)+": list(range("+str(numSpws)+"))"
            casaCmd = casaCmd + "}\n\n"

            casaCmd = casaCmd + "# spwmap for flux_inf and ampli_inf: matched by baseband to preserve signal path and keep science SPWs statistically independent as in non-BWSW case\n"
            casaCmd = casaCmd + "calspwmapf = "+repr(calspwmapf)+"\n\n"

            gainTableBase = gainTable.copy()
            gainTable.append(phaseForCal+str(bpassCalId))

            casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  field = '"+str(bpassCalId)+"', # "+fieldNames[bpassCalId]+"\n"
            casaCmd = casaCmd + "  gaintable = "+str(gainTable)+",\n"
            casaCmd = casaCmd + "  gainfield = ['"+str(bpassCalId)+"', '', '"+str(bpassCalId)+"'],\n"
            casaCmd = casaCmd + "  interp = [],\n"

            casaCmd = casaCmd + "  spwmap = [[], [], calspwmap["+str(bpassCalId)+"]],\n"

            casaCmd = casaCmd + "  calwt = True,\n"
            casaCmd = casaCmd + "  flagbackup = False)\n\n"

            gainTable.append(ampForSci[0])

            for i in calFieldIds:

                if i == bpassCalId: continue

                gainTable = gainTableBase.copy()
                gainTable.append(phaseForCal+str(i))
                gainTable.append(ampForSci[0])
                gainTable.append(flux)

                casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
                casaCmd = casaCmd + "  field = '"+str(i)+"', # "+fieldNames[i]+"\n"
                casaCmd = casaCmd + "  gaintable = "+str(gainTable)+",\n"

                casaCmd = casaCmd + "  gainfield = ['', '', '"+str(i)+"', '"+str(i)+"', '"+str(i)+"'],\n"
                casaCmd = casaCmd + "  interp = ['', 'nearest', 'linearPD', '', 'linear'],\n"
                casaCmd = casaCmd + "  spwmap = [[], [], calspwmap["+str(i)+"], calspwmapf, calspwmapf],\n"

                casaCmd = casaCmd + "  calwt = True,\n"
                casaCmd = casaCmd + "  flagbackup = False)\n\n"

    gainTable = []
    gainTable.append(bandpass)

    gainTableInterp = "'linear'"

    if len(phaseDiffCalTableName) != 0:

        if len(ampForSci) != 1: sys.exit('ERROR: missing table')

        gainTable.append(phaseDiffCalTableName[0])

        gainTable.append(phaseForSci+str(phaseCalFieldIds[0]))
        gainTable.append(ampForSci[0])
        gainTable.append(flux)
        gainTableInterp = "['', 'nearest', 'linearPD', '', 'linear']"

    else:

        gainTable.append(phaseForSci)
        gainTable.append(flux)
        gainTableInterp = "['linear','linear','linear']"
    
    for i in phaseCalFieldIds:

        sciFieldIds = []
        sciFieldNames = []
        for j in phaseCal:
            if phaseCal[j]['phaseCalId'] == i:
                sciFieldIds1 = phaseCal[j]['sciFieldIds']
                for k in sciFieldIds1: sciFieldIds.append(k)
                sciFieldNames.append(j)
        sciFieldIds = sorted(sciFieldIds)
        sciFieldNames = ','.join(sciFieldNames)

        if len(sciFieldIds) > 1:
            j0 = 0
            sciFieldIds1 = str(sciFieldIds[j0])
            for j in range(len(sciFieldIds)-1):
                if sciFieldIds[j+1] == sciFieldIds[j]+1: continue
                sciFieldIds1 = sciFieldIds1 + '~' + str(sciFieldIds[j])
                j0 = j+1
                sciFieldIds1 = sciFieldIds1 + ',' + str(sciFieldIds[j0])
            sciFieldIds1 = sciFieldIds1 + '~' + str(sciFieldIds[j+1])
        else:
            sciFieldIds1 = str(sciFieldIds[0])

        if renorm=='': # no renorm caltable
   
            casaCmd = casaCmd + "\napplycal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  field = '"+str(i)+","+sciFieldIds1+"', # "+sciFieldNames+"\n"
            casaCmd = casaCmd + "  gaintable = "+str(gainTable)+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "  gainfield = ['', '', '"+str(i)+"', '"+str(i)+"', '"+str(i)+"'], # "+fieldNames[i]+"\n"
            else:
                casaCmd = casaCmd + "  gainfield = ['', '"+str(i)+"', '"+str(i)+"'], # "+fieldNames[i]+"\n"
            casaCmd = casaCmd + "  interp = "+gainTableInterp+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "  spwmap = [[], [], calspwmap["+str(i)+"], calspwmapf, calspwmapf],\n"

            casaCmd = casaCmd + "  calwt = True,\n"
            casaCmd = casaCmd + "  flagbackup = False)\n"

        else: # we have a renormalization caltable; need to separate phasecalfields from scifields
            
            casaCmd = casaCmd + "\nif applyRenorm:\n"
            casaCmd = casaCmd + "  applycal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "    field = '"+str(i)+"', # phasecal\n"
            casaCmd = casaCmd + "    gaintable = "+str(gainTable)+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "    gainfield = ['', '', '"+str(i)+"', '"+str(i)+"', '"+str(i)+"'], # "+fieldNames[i]+"\n"
            else:
                casaCmd = casaCmd + "    gainfield = ['', '"+str(i)+"', '"+str(i)+"'], # "+fieldNames[i]+"\n"
            casaCmd = casaCmd + "    interp = "+gainTableInterp+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "    spwmap = [[], [], calspwmap["+str(i)+"], calspwmapf, calspwmapf],\n"

            casaCmd = casaCmd + "    calwt = True,\n"
            casaCmd = casaCmd + "    flagbackup = False)\n\n"

            gtable2 = gainTable.copy()
            gtable2.append(renorm)
            if len(phaseDiffCalTableName) != 0:
                gainTableInterp2 = "['', 'nearest', 'linearPD', '', 'linear', 'linear']"
            else:
                gainTableInterp2 = "['linear','linear','linear','linear']"

            casaCmd = casaCmd + "\n  applycal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "    field = '"+sciFieldIds1+"', # "+sciFieldNames+"\n"
            casaCmd = casaCmd + "    gaintable = "+str(gtable2)+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "    gainfield = ['', '', '"+str(i)+"', '"+str(i)+"', '"+str(i)+"', ''], # "+fieldNames[i]+"\n"
            else:
                casaCmd = casaCmd + "    gainfield = ['', '"+str(i)+"', '"+str(i)+"', ''], # "+fieldNames[i]+"\n"
            casaCmd = casaCmd + "    interp = "+gainTableInterp2+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "    spwmap = [[], [], calspwmap["+str(i)+"], calspwmapf, calspwmapf, []],\n"

            casaCmd = casaCmd + "    calwt = True,\n"
            casaCmd = casaCmd + "    flagbackup = False)\n"

            # also cover the case that renorm was found _not_ to be necessary
            
            casaCmd = casaCmd + "\nelse: # no renorm\n"
            casaCmd = casaCmd + "  applycal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "    field = '"+str(i)+","+sciFieldIds1+"', # "+sciFieldNames+"\n"
            casaCmd = casaCmd + "    gaintable = "+str(gainTable)+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "    gainfield = ['', '', '"+str(i)+"', '"+str(i)+"', '"+str(i)+"'], # "+fieldNames[i]+"\n"
            else:
                casaCmd = casaCmd + "    gainfield = ['', '"+str(i)+"', '"+str(i)+"'], # "+fieldNames[i]+"\n"
            casaCmd = casaCmd + "    interp = "+gainTableInterp+",\n"
            if len(phaseDiffCalTableName) != 0:
                casaCmd = casaCmd + "    spwmap = [[], [], calspwmap["+str(i)+"], calspwmapf, calspwmapf],\n"

            casaCmd = casaCmd + "    calwt = True,\n"
            casaCmd = casaCmd + "    flagbackup = False)\n"


            
    return casaCmd

#####################################
def doApplyB2BBandpassAndGainCalTables(msName, msName1='', iHaveSplitMyScienceSpw=False, bandpass='', renorm='', valueMaps={},
                                       combineB2BLFspws=False, combineB2BLFHFspws=False, combineB2BDGCspws=False):
    """Generate code for the final applycal step for the given MS (B2B case)"""

    print('\n*** doApplyB2BBandpassAndGainCalTables ***')
    print('Gathering information ...')

    if msName1 == '': msName1 = msName

    if bandpass == '': 
        casalog.post('ERROR: No bandpass cal table specified.', 'SEVERE')
        return False

    ### determine the scan lists and the field and spw ids to be used in the code

    if msName in valueMaps.keys():
        vm = valueMaps[msName]
        print('Using canned ValueMap.')
    else:
        vm = aU.ValueMapping(msName)
        valueMaps[msName] = vm

    intentSources = sfsdr.getIntentsAndSourceNames(msName)

    targetIds = intentSources['OBSERVE_TARGET']['id']
    targetIdStr = ','.join([str(i) for i in targetIds])
    checkSourceId  = intentSources['OBSERVE_CHECK']['id'][0]
    diffGainCalId = intentSources['CALIBRATE_DIFFGAIN']['id'][0]
    phaseCalId  = intentSources['CALIBRATE_PHASE']['id'][0]
    i = 0
    while phaseCalId == diffGainCalId and i < len(intentSources['CALIBRATE_PHASE']['id'])-1:
        i += 1
        phaseCalId  = intentSources['CALIBRATE_PHASE']['id'][i]

    fluxCalId  = intentSources['CALIBRATE_FLUX']['id'][0]

    if intentSources['CALIBRATE_FLUX']['spw'] == intentSources['CALIBRATE_DIFFGAIN']['spw']: # fluxcal covers HF and LF
        fluxCalLFId = fluxCalId
    else:
        fluxCalLFId = diffGainCalId 

    bpCalId  = intentSources['CALIBRATE_BANDPASS']['id'][0]

    spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#SIGNAL')
    if spwsDiffgainsig == []:
        spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#ON_SOURCE')
        sigintent = 'CALIBRATE_DIFFGAIN#ON_SOURCE'
    else:
        sigintent = 'CALIBRATE_DIFFGAIN#SIGNAL'

    spwHighDict = sfsdr.getSpwInfo(msName,intent=sigintent, caching=True)
    spwLowDict = sfsdr.getSpwInfo(msName,intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
    spwHigh = sorted(spwHighDict.keys())
    spwLow = sorted(spwLowDict.keys())

    mymsmd = msmdtool()
    mymsmd.open(msName)
    bpScans = mymsmd.scansforintent('CALIBRATE_BANDPASS#ON_SOURCE')
    bpCalScanList = ','.join([str(i) for i in bpScans])
    fluxCalScanList = ','.join([str(i) for i in mymsmd.scansforintent('CALIBRATE_FLUX#ON_SOURCE')])
    diffGainCalScanListLow  = ','.join([str(i) for i in mymsmd.scansforintent('CALIBRATE_DIFFGAIN#REFERENCE')])
    diffGainCalScanListHigh = ','.join([str(i) for i in mymsmd.scansforintent(sigintent)])
    lfScans = mymsmd.scansforspw(spwLow[0])
    mymsmd.close()


    mytb = aU.createCasaTool(tbtool)
    mytb.open(msName+'/OBSERVATION')
    obsTimeRange = mytb.getcol('TIME_RANGE')
    mytb.close()
    obsTimeStart = ((obsTimeRange[0]/86400.0)+2400000.5-2440587.5)*86400.0
    obsTimeStart = timeUtilities.strftime('%Y-%m-%dT%H:%M:%S', timeUtilities.gmtime(int(obsTimeStart[0])))



    # determine spw mapping 
    if iHaveSplitMyScienceSpw == True:
        spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
        spwIds = sorted(spwInfo.keys())

        spwHighSplit = []
        for myspw in spwHigh:
            spwHighSplit.append(spwIds.index(myspw))
        spwHighStr = ','.join([str(i) for i in spwHighSplit])

        spwLowSplit = []
        for myspw in spwLow:
            spwLowSplit.append(spwIds.index(myspw))
        spwLowStr = ','.join([str(i) for i in spwLowSplit])

        LFtoHF = list(range(max(spwHighSplit)+1))
        LFtoLF  = list(range(max(spwLowSplit)+1))
        HFtoHF  = list(range(max(spwHighSplit)+1))
        HFtoHF_DGC  = list(range(max(spwHighSplit)+1))
                    

        lowestSpwLow = int(min(spwLowSplit))
        
        if obsTimeStart > '2021-03-01T00:00:00':
            for spwRep in spwHighSplit: 
                BBnum = spwHighDict[spwRep]['basebandNum'] 
                for spwMap in spwLowSplit:  
                    if spwLowDict[spwMap]['basebandNum']==BBnum:
                        if combineB2BLFspws: 
                            LFtoHF[spwRep] = lowestSpwLow
                        else:
                            LFtoHF[spwRep] = spwMap
                            
        else:
            if combineB2BLFspws: 
                for spwRep in spwHighSplit:
                    LFtoHF[spwRep] = lowestSpwLow
            else:
                for spwId,spwRep in enumerate(spwHighSplit):
                    LFtoHF[spwRep] = spwLowSplit[spwId]

        # LFtoLF, HFtoHF, and HFtoHF_DGC are 'easy' mappings
        if combineB2BLFspws:
            for spwMap in spwLowSplit:
                LFtoLF[spwMap] = min(spwLowSplit)  # maps to min of HF index
        if combineB2BLFHFspws:
            for spwMap in spwHighSplit:
                HFtoHF[spwMap] = min(spwHighSplit)  # maps to min of HF index
        if combineB2BDGCspws:
            for spwMap in spwHighSplit:
                HFtoHF_DGC[spwMap] = min(spwHighSplit)  # maps to min of HF index

    else: # no reindexing took place
        spwHighStr = ','.join([str(i) for i in spwHigh])
        spwLowStr = ','.join([str(i) for i in spwLow])

        LFtoHF = list(range(max(spwHigh)+1))
        LFtoLF  = list(range(max(spwLow)+1))
        HFtoHF  = list(range(max(spwHigh)+1))
        HFtoHF_DGC  = list(range(max(spwHigh)+1))
        
        lowestSpwLow = int(min(spwLow))

        if obsTimeStart > '2021-03-01T00:00:00':
            for spwRep in spwHigh: 
                BBnum = spwHighDict[spwRep]['basebandNum'] 
                for spwMap in spwLow:  
                    if spwLowDict[spwMap]['basebandNum']==BBnum:
                        if combineB2BLFspws:
                             LFtoHF[spwRep] = lowestSpwLow
                        else:
                            LFtoHF[spwRep] = spwMap
        else:
            if combineB2BLFspws:
                for spwRep in spwHigh:
                    LFtoHF[spwRep] = lowestSpwLow
            else:
                for spwId,spwRep in enumerate(spwHigh):
                    LFtoHF[spwRep] = spwLow[spwId]

        # LFtoLF, HFtoHF, and HFtoHF_DGC are 'easy' mappings
        if combineB2BLFspws:
            for spwMap in spwLow:
                LFtoLF[spwMap] = min(spwLow)  # maps to min of HF index
        if combineB2BLFHFspws:
            for spwMap in spwHigh:
                HFtoHF[spwMap] = min(spwHigh)  # maps to min of HF index
        if combineB2BDGCspws:
            for spwMap in spwHigh:
                HFtoHF_DGC[spwMap] = min(spwHigh)  # maps to min of HF index


    print('bpCalScanList ', bpCalScanList)
    print('fluxCalScanList ', fluxCalScanList)
    print('diffGainCalScanListLow ', diffGainCalScanListLow)
    print('diffGainCalScanListHigh ', diffGainCalScanListHigh)
    print('spwHighStr ', spwHighStr)
    print('spwLowStr ', spwLowStr)

    print('LFtoHF ', LFtoHF)
    if combineB2BLFspws:
        print('LFtoLF ', LFtoLF)
    if combineB2BLFHFspws:
        print('HFtoHF ', HFtoHF)
    if combineB2BDGCspws:
        print('HFtoHF_DGC ', HFtoHF_DGC)
    ###

    print('Writing code ...')

    casaCmd = ''

    if bpCalId != fluxCalId:

        casaCmd = casaCmd + "# calibrating the bandpass calibrator\n"
        casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  field = '"+str(bpCalId)+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
        casaCmd = casaCmd + "  scan = '"+bpCalScanList+"',\n"
        casaCmd = casaCmd + "  interp=['linear','linear','linear','linear'],\n"
        if combineB2BLFHFspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(HFtoHF)+",[]],\n"
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phasehigh_int','"+msName1+".ampli_inf'],\n"
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"', '', '"+str(bpCalId)+"', '"+str(fluxCalId)+"'])\n"

        casaCmd = casaCmd + "\n# calibrating the fluxcalibrator\n"
        casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  field = '"+str(fluxCalId)+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n" 
        casaCmd = casaCmd + "  scan = '"+fluxCalScanList+"',\n"
        casaCmd = casaCmd + "  interp=['linear','linear','linear','linear'],\n"
        if combineB2BLFHFspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(HFtoHF)+",[]],\n"
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phasehigh_int','"+msName1+".ampli_inf'],\n"
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"', '','"+str(fluxCalId)+"', '"+str(fluxCalId)+"'])\n"

    else:
        casaCmd = casaCmd + "# calibrating the bandpass/fluxcalibrator\n"
        casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  field = '"+str(bpCalId)+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
        casaCmd = casaCmd + "  scan = '"+bpCalScanList+"',\n"
        casaCmd = casaCmd + "  interp=['linear','linear','linear''linear'],\n"
        if combineB2BLFHFspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(HFtoHF)+",[]],\n"        
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phasehigh_int','"+msName1+".ampli_inf'],\n"
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"', '', '"+str(bpCalId)+"', '"+str(fluxCalId)+"'])\n"

    # DGC is corrected with the LF phases
    # note we need to specify the scans, because if the DGC = BP, then the
    # previous BP correction is overwritten

    casaCmd = casaCmd + "\n# calibrating the DGC with the LF phases\n"
    casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
    casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
    casaCmd = casaCmd + "  scan = '"+diffGainCalScanListHigh+"',\n"
    casaCmd = casaCmd + "  interp=['linear','linear','linearPD','linear','linear'],\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_inf','"+msName1+".DGCoffset_multi','"+msName1+".ampli_inf'],\n"
    if combineB2BDGCspws:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoHF)+","+str(HFtoHF_DGC)+",[]],\n"
    else:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoHF)+",[],[]],\n"
    casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','','"+str(diffGainCalId)+"','"+str(diffGainCalId)+"','"+str(fluxCalId)+"'])\n"


    if renorm=='': # no renormalisation
        casaCmd = casaCmd + "\n## calibrating the target\n"
        casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  field = '"+targetIdStr+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n" 
        casaCmd = casaCmd + "  interp = ['linear','linear','linearPD', 'linear','linear'],\n"
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_inf','"+msName1+".DGCoffset_multi','"+msName1+".ampli_inf'],\n"
        if combineB2BDGCspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoHF)+","+str(HFtoHF_DGC)+",[]],\n"
        else:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoHF)+",[],[]],\n"        
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"', '', '"+str(phaseCalId)+"','"+str(diffGainCalId)+"','"+str(fluxCalId)+"'])\n"
    else: # with renorm caltable
        casaCmd = casaCmd + "\n## calibrating the target\n"
        casaCmd = casaCmd + "if applyRenorm:\n"
        casaCmd = casaCmd + "  applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "    field = '"+targetIdStr+"',\n"
        casaCmd = casaCmd + "    spw = '"+spwHighStr+"',\n" 
        casaCmd = casaCmd + "    interp = ['linear','linear','linearPD', 'linear','linear','linear'],\n"
        casaCmd = casaCmd + "    gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_inf','"+msName1+".DGCoffset_multi','"+msName1+".ampli_inf','"+renorm+"'],\n"
        if combineB2BDGCspws:
            casaCmd = casaCmd + "    spwmap = [[],[],"+str(LFtoHF)+","+str(HFtoHF_DGC)+",[],[]],\n"
        else:
            casaCmd = casaCmd + "    spwmap = [[],[],"+str(LFtoHF)+",[],[],[]],\n"
        casaCmd = casaCmd + "    gainfield = ['"+str(bpCalId)+"', '', '"+str(phaseCalId)+"','"+str(diffGainCalId)+"','"+str(fluxCalId)+"',''])\n"
        casaCmd = casaCmd + "else:\n"
        casaCmd = casaCmd + "  applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "    field = '"+targetIdStr+"',\n"
        casaCmd = casaCmd + "    spw = '"+spwHighStr+"',\n" 
        casaCmd = casaCmd + "    interp = ['linear','linear','linearPD', 'linear','linear'],\n"
        casaCmd = casaCmd + "    gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_inf','"+msName1+".DGCoffset_multi','"+msName1+".ampli_inf'],\n"
        if combineB2BDGCspws:
            casaCmd = casaCmd + "    spwmap = [[],[],"+str(LFtoHF)+","+str(HFtoHF_DGC)+",[]],\n"
        else:
            casaCmd = casaCmd + "    spwmap = [[],[],"+str(LFtoHF)+",[],[]],\n"
        casaCmd = casaCmd + "    gainfield = ['"+str(bpCalId)+"', '', '"+str(phaseCalId)+"','"+str(diffGainCalId)+"','"+str(fluxCalId)+"'])\n"

    casaCmd = casaCmd + "\n## calibrating the checksource\n"
    casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  field = '"+str(checkSourceId)+"',\n"
    casaCmd = casaCmd + "  spw = '"+spwHighStr+"',\n"
    casaCmd = casaCmd + "  interp = ['linear','linear','linearPD', 'linear','linear'],\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_inf','"+msName1+".DGCoffset_multi','"+msName1+".ampli_inf'],\n"
    if combineB2BDGCspws:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoHF)+","+str(HFtoHF_DGC)+",[]],\n"
    else:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoHF)+",[],[]],\n"        
    casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"', '', '"+str(phaseCalId)+"','"+str(diffGainCalId)+"','"+str(fluxCalId)+"'])\n"

    
    # LF data calibration

    casaCmd = casaCmd + "\n## calibration of the LF data\n"

    if bpCalId != fluxCalId and bpCalId != diffGainCalId and fluxCalId != diffGainCalId:

        casaCmd = casaCmd + "# calibrating the bandpass calibrator LF\n"
        casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  field = '"+str(bpCalId)+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
        casaCmd = casaCmd + "  scan = '"+bpCalScanList+"',\n"
        casaCmd = casaCmd + "  interp=['linear','linear','linear','linear'],\n"
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_int','"+msName1+".ampli_inf'],\n"
        if combineB2BLFspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoLF)+",[]],\n"
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','', '"+str(bpCalId)+"', '"+str(fluxCalId)+"'])\n"

        casaCmd = casaCmd + "\n# calibrating the fluxcalibrator LF\n"
        casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
        casaCmd = casaCmd + "  field = '"+str(fluxCalId)+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n" 
        casaCmd = casaCmd + "  scan = '"+fluxCalScanList+"',\n"
        casaCmd = casaCmd + "  interp=['linear','linear','linear','linear'],\n"
        casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_int','"+msName1+".ampli_inf'],\n"
        if combineB2BLFspws:
            casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoLF)+",[]],\n"
        casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','', '"+str(fluxCalId)+"', '"+str(fluxCalId)+"'])\n"

    else:
        #Test if the bandpass scans contain LF SPWs
        found=False
        for myscan in bpScans:
            if myscan in lfScans:
                found=True

        if found: # the BP scans do contain LF SPWs
            casaCmd = casaCmd + "# calibrating the bandpass/fluxcalibrator LF\n"
            casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
            casaCmd = casaCmd + "  field = '"+str(bpCalId)+"',\n"
            casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
            casaCmd = casaCmd + "  scan = '"+bpCalScanList+"',\n"
            casaCmd = casaCmd + "  interp=['linear','linear','linear','linear'],\n"
            casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_int','"+msName1+".ampli_inf'],\n"
            if combineB2BLFspws:
                casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoLF)+",[]],\n"
            casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','', '"+str(bpCalId)+"', '"+str(fluxCalId)+"'])\n"

    
    casaCmd = casaCmd + "\n# calibrating the DGC LF\n"
    casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  field = '"+str(diffGainCalId)+"',\n"
    casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
    casaCmd = casaCmd + "  scan = '"+diffGainCalScanListLow+"',\n"
    casaCmd = casaCmd + "  interp=['linear','linear','linear','linear'],\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_inf','"+msName1+".ampli_inf'],\n"
    if combineB2BLFspws:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoLF)+",[]],\n"
    casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','','"+str(diffGainCalId)+"','"+str(fluxCalLFId)+"'])\n"

    casaCmd = casaCmd + "\n# calibrating the phasecal LF\n"
    casaCmd = casaCmd + "applycal(vis = '"+msName1+"',\n"
    casaCmd = casaCmd + "  field = '"+str(phaseCalId)+"',\n"
    casaCmd = casaCmd + "  spw = '"+spwLowStr+"',\n"
    casaCmd = casaCmd + "  interp=['linear','linear','linear','linear'],\n"
    casaCmd = casaCmd + "  gaintable = ['"+bandpass+"','"+msName1+".DGCphase_phasediff','"+msName1+".phaselow_inf','"+msName1+".ampli_inf'],\n"
    if combineB2BLFspws:
        casaCmd = casaCmd + "  spwmap = [[],[],"+str(LFtoLF)+",[]],\n"
    casaCmd = casaCmd + "  gainfield = ['"+str(bpCalId)+"','','"+str(phaseCalId)+"','"+str(fluxCalLFId)+"'])\n"
    

    return casaCmd

#####################################
def doRenorm(msName, msName1='', isB2B=False, isBWSW=False, iHaveSplitMyScienceSpw=False, valueMaps={}):
    """Generate code for the renorm step needed for FDM data"""

    print('\n*** doRenorm ***')
    print('Gathering information ...')

    if msName1 == '': msName1 = msName

    spwRenorm = [] # default: let renorm code determine this automatically

    if isB2B or isBWSW: # revise spwRenorm
        print('Determining SPWs to renormalise ...')
        sciSpwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET', caching=True)
        spwInfo = sciSpwInfo.copy()
        #refSpwInfo = sfsdr.getSpwInfo(msName, intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
        #spwInfo.update(refSpwInfo)

        for myspw in sorted(spwInfo.keys()):
            if spwInfo[myspw]['numChans'] != 128: # must be FDM
                spwRenorm.append(myspw)

        if iHaveSplitMyScienceSpw == True:
            spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
            spwIds = sorted(spwInfo.keys())
            spwRenormSplit = []
            for myspw in spwRenorm:
                spwRenormSplit.append(spwIds.index(myspw))
            spwRenorm = spwRenormSplit

        print('spwRenorm: '+str(spwRenorm))

    ###

    print('Writing code ...')

    casaCmd = ''

    casaCmd = casaCmd + "# For background on this step see:\n"
    casaCmd = casaCmd + "# https://help.almascience.org/kb/articles/what-are-the-amplitude-calibration-issues-caused-by-alma-s-normalization-strategy\n\n"

    casaCmd = casaCmd + "from pipeline.extern.almarenorm import ACreNorm\n\n"

    casaCmd = casaCmd + "# QA2 analyst: First run this step individually with the following variable *applyRenorm* set to False\n"
    casaCmd = casaCmd + "#              in order to determine the max. renormalization scaling factor.\n"
    casaCmd = casaCmd + "#              If the factor is <= 1.02, don't change anything and proceed with the remaining steps.\n"
    casaCmd = casaCmd + "#              If the max. renormalization scaling factor is > 1.02, then edit the following line in this script\n"
    casaCmd = casaCmd + "#              and set applyRenorm=True, and then re-run this step.\n"
    casaCmd = casaCmd + "#              In any case, note here the max. renorm factor you determined:  <max. renorm factor determined in this step>\n"
    casaCmd = casaCmd + "applyRenorm = False\n\n"

    casaCmd = casaCmd + "if applyRenorm or applyonly!=True:\n\n"

    casaCmd = casaCmd + "  if not applyRenorm: os.system('rm -rf RN_plots')\n\n" 
    casaCmd = casaCmd + "  rn=ACreNorm('"+msName1+"')\n\n"
    casaCmd = casaCmd + "  corrApplied = rn.checkApply()\n" 
    casaCmd = casaCmd + "  if corrApplied:\n"
    casaCmd = casaCmd + "    casalog.post(' The renormalization has already been applied to these data', 'WARN')\n"
    casaCmd = casaCmd + "    if applyRenorm:\n"
    casaCmd = casaCmd + "        sys.exit('Do not rerun the renorm step again. Renorm was applied already.')\n\n"

    casaCmd = casaCmd + "  if applyRenorm or applyonly==True:\n"
    casaCmd = casaCmd + "    diagSpectra=False\n"
    casaCmd = casaCmd + "  else:\n"
    casaCmd = casaCmd + "    diagSpectra=True\n\n"
    
    if aU.getCasaVersion() < "6.4.1":
        adtlRenormSettings = ""
    else:
        adtlRenormSettings = ", atmAutoExclude=True"

    if isB2B or isBWSW:
        casaCmd = casaCmd + "  rn.renormalize(docorr=applyRenorm, spws="+str(spwRenorm)+", diagSpectra=diagSpectra, antHeuristicsSpectra=False"+adtlRenormSettings+")\n\n"
    else:
        casaCmd = casaCmd + "  rn.renormalize(docorr=applyRenorm, diagSpectra=diagSpectra, antHeuristicsSpectra=False"+adtlRenormSettings+")\n\n"

    casaCmd = casaCmd + "  if applyRenorm or applyonly==True:\n"
    casaCmd = casaCmd + "    pass  # no plotting\n"
    casaCmd = casaCmd + "  else:\n"
    casaCmd = casaCmd + "    rn.plotSpectra()\n\n"

    casaCmd = casaCmd + "  if os.path.exists('RN_plots') and applyRenorm:\n" 
    casaCmd = casaCmd + "    os.system('rm -f RN_plots/*.pdf RN_plots/*_scan*_field*.png')\n" 
    casaCmd = casaCmd + "    os.system('mv RN_plots "+msName1+".renorm.plots')\n\n" 

    casaCmd = casaCmd + "  if not applyRenorm:\n" 

    if aU.getCasaVersion() < "6.4.1":
        casaCmd = casaCmd + "    maxScaling = rn.MaxOutScaling\n"
    else:
        casaCmd = casaCmd + "    mystats = rn.rnpipestats\n"
        casaCmd = casaCmd + "    casalog.post('Maximum Renormalization scaling factors:\\n'+str(mystats))\n\n"  
        casaCmd = casaCmd + "    rnfactors = []\n" 
        casaCmd = casaCmd + "    for field in mystats.keys():\n" 
        casaCmd = casaCmd + "      print(' ************************')\n"
        casaCmd = casaCmd + "      print(field, ': max. renorm. scaling for each SPW')\n" 
        casaCmd = casaCmd + "      for spw in mystats[field].keys():\n" 
        casaCmd = casaCmd + "        print(spw,': ', mystats[field][spw])\n"     
        casaCmd = casaCmd + "        rnfactors.append(mystats[field][spw]['max_rn'])\n"     
        casaCmd = casaCmd + "    maxScaling = max(rnfactors)\n"

    casaCmd = casaCmd + "    print(' ************************')\n"
    casaCmd = casaCmd + "    print(' The maximum Renormalization scaling is '+str(maxScaling))\n"
    casaCmd = casaCmd + "    if maxScaling>1.02:\n"
    casaCmd = casaCmd + "      casalog.post('Renormalization should be applied before proceeding.','WARN')\n"
    casaCmd = casaCmd + "      rn.close()\n"
    casaCmd = casaCmd + "      sys.exit('Please edit the renorm step to set applyRenorm=True and rerun the step!')\n"
    casaCmd = casaCmd + "    else:\n"
    casaCmd = casaCmd + "      print(' No Renormalization needed.')\n"
    casaCmd = casaCmd + "    print(' ************************')\n\n"

    casaCmd = casaCmd + "  rn.close()\n"

    return casaCmd

#####################################
def doRenormTable(msName, msName1='', isB2B=False, isBWSW=False, iHaveSplitMyScienceSpw=False, valueMaps={}):
    """Generate code for the renorm step needed for FDM data"""

    print('\n*** doRenormTable ***')
    print('Gathering information ...')

    if msName1 == '': msName1 = msName
    renormTabName = msName1+'.renorm.tbl'

    spwRenorm = [] # default: let renorm code determine this automatically
    

    if isB2B or isBWSW: # revise spwRenorm
        print('Determining SPWs to renormalise ...')
        sciSpwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET', caching=True)
        spwInfo = sciSpwInfo.copy()
        #refSpwInfo = sfsdr.getSpwInfo(msName, intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)
        #spwInfo.update(refSpwInfo)

        for myspw in sorted(spwInfo.keys()):
            if spwInfo[myspw]['numChans'] != 128: # must be FDM
                spwRenorm.append(myspw)

        if iHaveSplitMyScienceSpw == True:
            spwInfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET|CALIBRATE_BANDPASS', caching=True)
            spwIds = sorted(spwInfo.keys())
            spwRenormSplit = []
            for myspw in spwRenorm:
                spwRenormSplit.append(spwIds.index(myspw))
            spwRenorm = spwRenormSplit

        print('spwRenorm: '+str(spwRenorm))

    ###

    print('Writing code ...')

    casaCmd = ''

    casaCmd = casaCmd + "# For background on this step see:\n"
    casaCmd = casaCmd + "# https://help.almascience.org/kb/articles/what-are-the-amplitude-calibration-issues-caused-by-alma-s-normalization-strategy\n\n"

    casaCmd = casaCmd + "os.system('rm -rf "+renormTabName+"')\n\n"

    casaCmd = casaCmd + "from pipeline.extern.almarenorm import ACreNorm\n\n"

    casaCmd = casaCmd + "# QA2 analyst: First run this step individually with the following variable *applyRenorm* set to False\n"
    casaCmd = casaCmd + "#              in order to determine the max. renormalization scaling factor.\n"
    casaCmd = casaCmd + "#              If the factor is <= 1.02, don't change anything and proceed with the remaining steps.\n"
    casaCmd = casaCmd + "#              If the max. renormalization scaling factor is > 1.02, then edit the following line in this script\n"
    casaCmd = casaCmd + "#              and set applyRenorm=True, and then re-run this step.\n"
    casaCmd = casaCmd + "#              In any case, note here the max. renorm factor you determined:  <max. renorm factor determined in this step>\n"
    casaCmd = casaCmd + "applyRenorm = False\n\n"

    casaCmd = casaCmd + "if applyRenorm or applyonly!=True:\n\n"

    casaCmd = casaCmd + "  if not applyRenorm: os.system('rm -rf RN_plots')\n\n" 
    casaCmd = casaCmd + "  rn=ACreNorm('"+msName1+"')\n\n"

    casaCmd = casaCmd + "  if applyRenorm or applyonly==True:\n"
    casaCmd = casaCmd + "    diagSpectra=False\n"
    casaCmd = casaCmd + "  else:\n"
    casaCmd = casaCmd + "    diagSpectra=True\n\n"
    
    if isB2B or isBWSW:
        casaCmd = casaCmd + "  rn.renormalize(createCalTable=applyRenorm, spws="+str(spwRenorm)+", diagSpectra=diagSpectra, antHeuristicsSpectra=False, atmAutoExclude=True, fillTable=True)\n\n"
    else:
        casaCmd = casaCmd + "  rn.renormalize(createCalTable=applyRenorm, diagSpectra=diagSpectra, antHeuristicsSpectra=False, atmAutoExclude=True, fillTable=True)\n\n"

    casaCmd = casaCmd + "  if applyRenorm or applyonly==True:\n"
    casaCmd = casaCmd + "    pass  # no plotting\n"
    casaCmd = casaCmd + "  else:\n"
    casaCmd = casaCmd + "    rn.plotSpectra()\n\n"

    casaCmd = casaCmd + "  if os.path.exists('RN_plots') and applyRenorm:\n" 
    casaCmd = casaCmd + "    os.system('rm -f RN_plots/*.pdf RN_plots/*_scan*_field*.png')\n" 
    casaCmd = casaCmd + "    os.system('mv RN_plots "+msName1+".renorm.plots')\n\n" 

    casaCmd = casaCmd + "  if not applyRenorm:\n" 

    casaCmd = casaCmd + "    mystats = rn.rnpipestats\n"
    casaCmd = casaCmd + "    casalog.post('Maximum Renormalization scaling factors:\\n'+str(mystats))\n\n"  
    casaCmd = casaCmd + "    rnfactors = []\n" 
    casaCmd = casaCmd + "    for field in mystats.keys():\n" 
    casaCmd = casaCmd + "      print(' ************************')\n"
    casaCmd = casaCmd + "      print(field, ': max. renorm. scaling for each SPW')\n" 
    casaCmd = casaCmd + "      for spw in mystats[field].keys():\n" 
    casaCmd = casaCmd + "        print(spw,': ', mystats[field][spw])\n"     
    casaCmd = casaCmd + "        rnfactors.append(mystats[field][spw]['max_rn'])\n"     
    casaCmd = casaCmd + "    maxScaling = max(rnfactors)\n"

    casaCmd = casaCmd + "    print(' ************************')\n"
    casaCmd = casaCmd + "    print(' The maximum Renormalization scaling is '+str(maxScaling))\n"
    casaCmd = casaCmd + "    if maxScaling>1.02:\n"
    casaCmd = casaCmd + "      casalog.post('Renormalization should be applied before proceeding.','WARN')\n"
    casaCmd = casaCmd + "      rn.close()\n"
    casaCmd = casaCmd + "      sys.exit('Please edit the renorm step to set applyRenorm=True and rerun the step!')\n"
    casaCmd = casaCmd + "    else:\n"
    casaCmd = casaCmd + "      print(' No Renormalization needed.')\n"
    casaCmd = casaCmd + "    print(' ************************')\n\n"

    casaCmd = casaCmd + "  rn.close()\n"

    return casaCmd, renormTabName


#####################################

def doFluxCalibration(msNames, fluxFile='allFluxes.txt', refant='', valueMaps={}):
    """Generate code for the 'flux equalisation' of a set of calibrated MSs."""

    print('\n*** doFluxCalibration ***')
    print('Gathering information ...')

    if type(msNames).__name__ == 'str': msNames = [msNames]

    if os.path.exists(fluxFile) == False: 
        casalog.post('ERROR: Flux file '+fluxFile+' does not seem to exist in the current directory.','SEVERE')
        return False

    casaCmd = 'print "# Flux calibration of the data."\n\n'

    f = open(fluxFile, 'r')
    fc = f.readlines()
    f.close()

    msName = []
    fieldName = []
    spwId = []
    fluxVal = []

    for line in fc:

        if len(line) == 0 or line.isspace() == True or line.lstrip()[0] == '#': continue

        casaCmd = casaCmd + '# ' + line

        line = line.split('"')
        fieldName.append(line[1])
        line = line[2].split()
        spwId.append(line[0])
        fluxVal.append(line[3])
        msName.append(line[5])

    msName = np.array(msName)
    fieldName = np.array(fieldName)

    casaCmd = casaCmd + '\n'

    for i in range(len(msName)):

        if msName[i] not in msNames: 
            casalog.post('ERROR: Missing dataset.','SEVERE')
            return False

        casaCmd = casaCmd + "setjy(vis = '"+msName[i]+"',\n"
        casaCmd = casaCmd + "  field = '"+fieldName[i]+"',\n"
        casaCmd = casaCmd + "  spw = '"+spwId[i]+"',\n"
        casaCmd = casaCmd + "  standard = 'manual',\n"
        casaCmd = casaCmd + "  fluxdensity = ["+fluxVal[i]+", 0, 0, 0])\n\n"

    for i in range(len(msNames)):

        calFieldNames = np.unique(fieldName[np.where(msName == msNames[i])])
        calFieldNames = ','.join(calFieldNames)

        myRefAnt = refant
        if myRefAnt == '': myRefAnt = sfsdr.getRefAntenna(msNames[i])
        casaCmd = casaCmd + "os.system('rm -rf %s.ampli_inf') \n"%(msNames[i])
        casaCmd = casaCmd + "gaincal(vis = '"+msNames[i]+"',\n"
        casaCmd = casaCmd + "  caltable = '"+msNames[i]+".ampli_inf',\n"
        casaCmd = casaCmd + "  field = '"+calFieldNames+"',\n"
        casaCmd = casaCmd + "  solint = 'inf',\n"
        casaCmd = casaCmd + "  combine = 'scan',\n"
        casaCmd = casaCmd + "  refant = '"+myRefAnt+"',\n"
        casaCmd = casaCmd + "  gaintype = 'T',\n"
        casaCmd = casaCmd + "  calmode = 'a')\n\n"

        phaseCal = sfsdr.getPhaseCal(msNames[i], valueMaps=valueMaps)

        for sciFieldName in list(phaseCal.keys()):

            sciFieldIds = phaseCal[sciFieldName]['sciFieldIds']

            if len(sciFieldIds) > 1:
                j0 = 0
                sciFieldIds1 = str(sciFieldIds[j0])
                for j in range(len(sciFieldIds)-1):
                    if sciFieldIds[j+1] == sciFieldIds[j]+1: continue
                    sciFieldIds1 = sciFieldIds1 + '~' + str(sciFieldIds[j])
                    j0 = j+1
                    sciFieldIds1 = sciFieldIds1 + ',' + str(sciFieldIds[j0])
                sciFieldIds1 = sciFieldIds1 + '~' + str(sciFieldIds[j+1])
            else:
                sciFieldIds1 = str(sciFieldIds[0])

            casaCmd = casaCmd + "applycal(vis = '"+msNames[i]+"',\n"
            casaCmd = casaCmd + "  field = '"+str(phaseCal[sciFieldName]['phaseCalId'])+","+sciFieldIds1+"', # "+phaseCal[sciFieldName]['phaseCalName']+","+sciFieldName+"\n"
            casaCmd = casaCmd + "  gaintable = '"+msNames[i]+".ampli_inf',\n"
            casaCmd = casaCmd + "  gainfield = '"+str(phaseCal[sciFieldName]['phaseCalId'])+"', # "+phaseCal[sciFieldName]['phaseCalName']+"\n"
            casaCmd = casaCmd + "  calwt = False,\n"
            casaCmd = casaCmd + "  flagbackup = False)\n\n"

    if len(msNames) > 1:
        casaCmd = casaCmd + 'print "# Concatenating the data."\n\n'
        casaCmd = casaCmd + "concat(vis = "+str([i for i in msNames])+",\n"
        casaCmd = casaCmd + "  concatvis = 'calibrated.ms')\n\n"

    return casaCmd

###################################

def SDdoFillTsysSolutions(asapName, msName='', spwIds='', tsysCalTableName='', tsysmap='', iHaveSplitMyScienceSpw=False, doplot=False, sky=False, calmode=''):
    """Generate code for the Tsys solution step of an SD calibration script.

    spwIds must be specified as a string, e.g. spwIds = '1,3,5,7'"""

    if msName == '':
        if spwIds == '' and tsysmap == '': 
            casalog.post('ERROR: you have not specified neither msName, or spwIds and tsysmap.','SEVERE')
            return False
        if doplot == True: 
            casalog.post('ERROR: you have not specified msName, so I cannot do any plot.', 'SEVERE')
            return False

    if type(asapName).__name__ == 'str': asapName = [asapName]

    casaCmd = ''

    if sky == False:

        calTableName1 = msName + '.tsys'

        casaCmd = casaCmd + "os.system('rm -Rf "+calTableName1+"')\n\n"
        casaCmd = casaCmd + "gencal(vis = '"+msName+"',\n"
        casaCmd = casaCmd + "  caltable = '"+calTableName1+"',\n"
        casaCmd = casaCmd + "  caltype = 'tsys')\n\n"

    else:

        calTableName1 = msName + '.sky'

        casaCmd = casaCmd + "os.system('rm -Rf "+calTableName1+"')\n\n"
        casaCmd = casaCmd + "sdcal(infile = '"+msName+"',\n"
        casaCmd = casaCmd + "  outfile = '"+calTableName1+"',\n"
        casaCmd = casaCmd + "  calmode = '"+calmode+"')\n\n"

    tsysCalTableName.append(calTableName1)

    if doplot == True:
        chanrange = '92.1875%'

        casaCmd = casaCmd + "plotbandpass(caltable='%s', overlay='time', \n" %(calTableName1)
        casaCmd = casaCmd + "  xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,\n"
        casaCmd = casaCmd + "  showatm=True,pwv='auto',chanrange='"+chanrange+"',showfdm=True, \n"
        casaCmd = casaCmd + "  field='', figfile='%s') \n" %(calTableName1+'.plots.overlayTime/'+calTableName1.split('/')[-1])

        casaCmd = casaCmd + "\nif applyonly != True: es.checkCalTable('"+calTableName1+"', msName='"+msName+"', interactive=False)\n"

    return casaCmd

#######################################

def SDdoCalibration(asapName, msName='', spwIds='', calmode='ps', tsysCalTableName='', tsysmap='', iHaveSplitMyScienceSpw=False, doplot=True, skyCalTableName='',
                    jyCalTableName=''):
    """Generate code for the sdcal step of an SD  calibration script.

       spwIds must be specified as a string, e.g. spwIds = '1,3,5,7'"""

    myCasaVersion = aU.getCasaVersion()

    if msName == '' and spwIds == '': 
        casalog.post('ERROR: you have not specified neither msName, or spwIds.','SEVERE')
        return False

    if type(asapName).__name__ == 'str': asapName = [asapName]

    if msName != '':
        spwInfo = sfsdr.getSpwInfo(msName, caching=True)
        spwIds1 = sorted(spwInfo.keys())
        spwIds1 = [int(i) for i in spwIds1]
        if iHaveSplitMyScienceSpw == True: 
            spwIds1 = list(range(len(spwIds1)))
        spwIds = ','.join([str(i) for i in spwIds1])
    else:
        spwIds1 = spwIds.split(',')
        spwIds1 = [int(i) for i in spwIds1]

    casaCmd = ''

    if tsysmap == '':
        if tsysCalTableName == '': 
            casalog.post('ERROR: you have not specified a Tsys cal table.','SEVERE')
            return False
        if myCasaVersion < '5.9.9':
            casaCmd = casaCmd + "from recipes.almahelpers import tsysspwmap\n"
        else:
            casaCmd = casaCmd + "from casarecipes.almahelpers import tsysspwmap\n"
        casaCmd = casaCmd + "tsysmap = tsysspwmap(vis = '"+msName+"', tsystable = '"+tsysCalTableName+"', trim = False)\n\n"
    else:
        casaCmd = casaCmd + "tsysmap = "+str(tsysmap)+"\n\n"

    if myCasaVersion < '5.0':
        mymsmd = msmdtool()
        mymsmd.open(msName)
        tsysspw = [i for i in mymsmd.spwsforintent('CALIBRATE_ATMOSPHERE#ON_SOURCE') if i not in mymsmd.chanavgspws().tolist()+mymsmd.wvrspws().tolist()]
        if (tsysspw == []):
            print("Will use CALIBRATE_ATMOSPHERE#HOT instead.")
            tsysspw = [i for i in mymsmd.spwsforintent('CALIBRATE_ATMOSPHERE#HOT') if i not in mymsmd.chanavgspws().tolist()+mymsmd.wvrspws().tolist()]

        mymsmd.close()

        casaCmd = casaCmd + "spwmap = {}\n"
        casaCmd = casaCmd + "for i in "+str(spwIds1)+":\n"
        casaCmd = casaCmd + "  if not tsysmap[i] in spwmap.keys():\n"
        casaCmd = casaCmd + "    spwmap[tsysmap[i]] = []\n"
        casaCmd = casaCmd + "  spwmap[tsysmap[i]].append(i)\n\n"

    for i in asapName:

        if myCasaVersion < '5.0':

            casaCmd = casaCmd + "os.system('rm -Rf "+i+".cal')\n\n"

            calmode1 =  calmode+',tsys,apply'

            casaCmd = casaCmd + "sdcal2(infile = '"+i+"',\n"
            casaCmd = casaCmd + "  calmode = '"+calmode1+"',\n"
            casaCmd = casaCmd + "  spw = '"+','.join([str(j) for j in sorted(spwIds1+tsysspw)])+"',\n"
            casaCmd = casaCmd + "  tsysspw = '"+','.join([str(j) for j in tsysspw])+"',\n"
            casaCmd = casaCmd + "  spwmap = spwmap,\n"
            casaCmd = casaCmd + "  outfile = '"+i+".cal',\n"
            casaCmd = casaCmd + "  overwrite = True)\n\n"

        else:
            mymsmd = msmdtool()
            mymsmd.open(msName)
            targetFieldIds = mymsmd.fieldsforintent('OBSERVE_TARGET#ON_SOURCE')
            offFieldIds = mymsmd.fieldsforintent('OBSERVE_TARGET#OFF_SOURCE')
            fieldNames = mymsmd.fieldnames()
            mymsmd.close()

            mygainfield = 'str(i)' # use the targetfield itself by default
            if len(offFieldIds)>0:
                mygainfield = 'offfields[i]'

                # determine array of OFF fields
                mygainfields = list(range(0, max(targetFieldIds)+1))
                for j in targetFieldIds:
                    tfieldname = fieldNames[j]
                    for k in offFieldIds:
                        ofieldname = fieldNames[k]
                        if ofieldname.find('_OFF_') >= 0:
                            if ofieldname[:ofieldname.find('_OFF_')] == tfieldname:
                                print('Using field '+ofieldname+' as OFF-SOURCE field for '+tfieldname+'.')
                                mygainfields[j] = k
                                break
                    if mygainfields[j]==j:
                        casalog.post('Target field '+tfieldname+' does not seem to have an OFF-SOURCE field. Will use the target field itself as gainfield.', 'WARN')
                            
                casaCmd = casaCmd + "offfields = "+str([str(mygainfields[j]) for j in range(0, max(targetFieldIds)+1)])+"\n\n"

            casaCmd = casaCmd + "for i in "+str([j for j in targetFieldIds])+":\n"
            casaCmd = casaCmd + "  applycal(vis = '"+i+"',\n"
            casaCmd = casaCmd + "    applymode = 'calflagstrict',\n"
            casaCmd = casaCmd + "    spw = '"+','.join([str(j) for j in sorted(spwIds1)])+"',\n"
            casaCmd = casaCmd + "    field = str(i),\n"
            if len(jyCalTableName) > 0:
                casaCmd = casaCmd + "    gaintable = ['"+tsysCalTableName+"', '"+skyCalTableName+"', '"+jyCalTableName+"'],\n"
            else:
                casaCmd = casaCmd + "    gaintable = ['"+tsysCalTableName+"', '"+skyCalTableName+"'],\n"

            casaCmd = casaCmd + "    gainfield = ['nearest', "+mygainfield+"],\n"
            casaCmd = casaCmd + "    spwmap = tsysmap)\n\n"

        if doplot == True:
            if myCasaVersion < '5.0':
                casaCmd = casaCmd + "if applyonly != True: es.SDcheckSpectra(asapName='"+i+".cal', spwIds='"+spwIds+"', interactive=False)\n\n"
            else:
                casaCmd = casaCmd + "if applyonly != True: es.SDcheckSpectra(msName='"+i+"', spwIds='"+spwIds+"', intent='OBSERVE_TARGET#ON_SOURCE', interactive=False)\n\n"

    return casaCmd


##################################################

def SDdoAtmCor(msName, jyCalTableName):
    """Generate code for the sdatmcor step of an SD calibration script.
    """

    myCasaVersion = aU.getCasaVersion()

    if myCasaVersion < '6.4.4':
        casalog.post('sdatmcor only available for CASA >= 6.4.4.', 'SEVERE')
        return False

    if msName == '': 
        casalog.post('ERROR: you have not specified  msName','SEVERE')
        return False

    spwInfo = sfsdr.getSpwInfo(msName, caching=True)

    mymsmd = msmdtool()
    mymsmd.open(msName)
    targetFieldIds = mymsmd.fieldsforintent('OBSERVE_TARGET#ON_SOURCE')
    mymsmd.close()

    casaCmd = ''

    casaCmd = casaCmd + "sdatmcor(infile = '"+msName+"',\n"
    casaCmd = casaCmd + "         datacolumn = 'corrected',\n"
    casaCmd = casaCmd + "         outfile = '"+msName+".atmcor.atmtype1',\n"
    casaCmd = casaCmd + "         overwrite = True,\n"
    casaCmd = casaCmd + "         field = '"+','.join([str(i) for i in targetFieldIds])+"',\n"
    casaCmd = casaCmd + "         intent = 'OBSERVE_TARGET#ON_SOURCE',\n"
    casaCmd = casaCmd + "         outputspw = '"+','.join([str(i) for i in spwInfo.keys()])+"',\n"
    casaCmd = casaCmd + "         gainfactor = '"+jyCalTableName+"',\n"
    casaCmd = casaCmd + "         atmtype=1)\n"
    
    return casaCmd


##################################################

def SDdoBaselineSubtraction(asapName, msName='', spwIds='', iHaveSplitMyScienceSpw=False, doplot=True):
    """Generate code for the baseline subtraction step of an SD calibration script.

       spwIds must be specified as a string, e.g. spwIds = '1,3,5,7'"""

    if msName == '' and spwIds == '': 
        casalog.post('ERROR: you have not specified neither msName, or spwIds.','SEVERE')

    if type(asapName).__name__ == 'str': asapName = [asapName]

    if msName != '':
        spwInfo = sfsdr.getSpwInfo(msName, caching=True)
        spwIds1 = sorted(spwInfo.keys())
        spwIds1 = [int(i) for i in spwIds1]
        if iHaveSplitMyScienceSpw == True: 
            spwIds1 = list(range(len(spwIds1)))
        spwIds = ','.join([str(i) for i in spwIds1])
    else:
        spwIds1 = spwIds.split(',')
        spwIds1 = [int(i) for i in spwIds1]

    casaCmd = ''

    for i in asapName:

        casaCmd = casaCmd + "os.system('rm -Rf "+i+".bl')\n\n"

        if aU.getCasaVersion() >= '5.0':

            if aU.getCasaVersion() < '6.4.4':

                casaCmd = casaCmd + "sdbaseline(infile = '"+i+"',\n"
                casaCmd = casaCmd + "  datacolumn = 'corrected',\n"

            else: # atmcor was applied beforehand

                casaCmd = casaCmd + "sdbaseline(infile = '"+i+".atmcor.atmtype1',\n"
                casaCmd = casaCmd + "  datacolumn = 'data',\n"

            casaCmd = casaCmd + "  spw = '"+','.join([str(j) for j in spwIds1])+"',\n"
            casaCmd = casaCmd + "  maskmode = 'auto',\n"
            casaCmd = casaCmd + "  thresh = 5.0,\n"
            casaCmd = casaCmd + "  avg_limit = 4,\n"
            casaCmd = casaCmd + "  blfunc = 'poly',\n"
            casaCmd = casaCmd + "  order = 1,\n"
            casaCmd = casaCmd + "  outfile = '"+i+".bl')\n\n"

        else:

            casaCmd = casaCmd + "sdbaseline(infile = '"+i+"',\n"
            casaCmd = casaCmd + "  spw = '"+','.join([str(j) for j in spwIds1])+"',\n"
            casaCmd = casaCmd + "  maskmode = 'auto',\n"
            casaCmd = casaCmd + "  thresh = 5.0,\n"
            casaCmd = casaCmd + "  avg_limit = 4,\n"
            casaCmd = casaCmd + "  blfunc = 'poly',\n"
            casaCmd = casaCmd + "  order = 1,\n"
            casaCmd = casaCmd + "  outfile = '"+i+".bl',\n"
            casaCmd = casaCmd + "  overwrite = True)\n\n"

        if doplot == True:
            if aU.getCasaVersion() < '5.0':
                casaCmd = casaCmd + "if applyonly != True: es.SDcheckSpectra(asapName='"+i+".bl', spwIds='"+spwIds+"', interactive=False)\n\n"
            else:
                spwIds = list(range(len(spwIds.split(','))))
                spwIds = ','.join([str(j) for j in spwIds])
                casaCmd = casaCmd + "if applyonly != True: es.SDcheckSpectra(msName='"+i+".bl', spwIds='"+spwIds+"', intent='OBSERVE_TARGET#ON_SOURCE', interactive=False)\n\n"

    return casaCmd

def isFullPol(msName, valueMaps={}):
    """Determine if the given MS is a full polarisation dataset.
    Return True if so."""

    if msName in valueMaps.keys():
        vm = valueMaps[msName]
        print('Using canned ValueMap.')
    else:
        vm = aU.ValueMapping(msName)
        valueMaps[msName] = vm

    polcalSpws = vm.getSpwsForIntent('CALIBRATE_POLARIZATION#ON_SOURCE')
    if polcalSpws == []:
        return False
    else:
        scienceSpws = vm.getSpwsForIntent('OBSERVE_TARGET#ON_SOURCE')
        for myspw in scienceSpws:
            if not myspw in polcalSpws:
                casalog.post("This dataset seems to be FULL POLARISATION (it has polcal data) but the polcal was not observed for all science SPWs.",
                             'WARN')
                break

    return True



def isB2BorBWSW(msName, valueMaps={}):
    """Determine if the given dataset uses band-to-band phase transfer
       or bandwidth switching. Returns two bool values
       isB2B, isBWSW"""

    if msName in valueMaps.keys():
        vm = valueMaps[msName]
        print('Using canned ValueMap.')
    else:
        vm = aU.ValueMapping(msName)
        valueMaps[msName] = vm

    spwsPhasecal = vm.getSpwsForIntent('CALIBRATE_PHASE#ON_SOURCE')
    spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#SIGNAL')
    if spwsDiffgainsig == []:
        spwsDiffgainsig = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#ON_SOURCE')
        diffGainSignalIntent = 'CALIBRATE_DIFFGAIN#ON_SOURCE'
    else:
        diffGainSignalIntent = 'CALIBRATE_DIFFGAIN#SIGNAL'
    spwsDiffgainref = vm.getSpwsForIntent('CALIBRATE_DIFFGAIN#REFERENCE')

    isB2B = False
    isBWSW = False

    if spwsPhasecal != []:
        print('spwsPhasecal '+str(spwsPhasecal))
        print('spwsDiffgainsig '+str(spwsDiffgainsig))
        print('spwsDiffgainref '+str(spwsDiffgainref))

        if (spwsDiffgainsig != []) and (spwsDiffgainref != []):
            print('Found intents '+diffGainSignalIntent+' and CALIBRATE_DIFFGAIN#REFERENCE.')
            casalog.post('Found intents '+diffGainSignalIntent+' and CALIBRATE_DIFFGAIN#REFERENCE.', 'INFO')

            if spwsPhasecal == spwsDiffgainref:
                minFreq=1E99
                maxFreq=0
                sigMaxbw=0
                refMaxbw=0
                for myspw in spwsDiffgainsig:
                    if myspw in vm.spwInfo.keys():
                        thefreq = vm.spwInfo[myspw]['meanFreq']
                        thebw = vm.spwInfo[myspw]['bandwidth']
                        thenumchan = vm.spwInfo[myspw]['numChannels']
                        if thefreq < minFreq:
                            minFreq = thefreq
                        if thefreq > maxFreq:
                            maxFreq = thefreq
                        if thebw > sigMaxbw and thenumchan > 4:
                            sigMaxbw = thebw

                for myspw in spwsDiffgainref:
                    if myspw in vm.spwInfo.keys():
                        thefreq = vm.spwInfo[myspw]['meanFreq']
                        thebw = vm.spwInfo[myspw]['bandwidth']
                        thenumchan = vm.spwInfo[myspw]['numChannels']
                        if thefreq < minFreq:
                            minFreq = thefreq
                        if thefreq > maxFreq:
                            maxFreq = thefreq
                        if thebw > refMaxbw and thenumchan > 4:
                            #if 'CALIBRATE_POINTING#ON_SOURCE' in vm.getIntentsForSpw(myspw):
                            #    print('    SPW '+str(myspw)+' has intent CALIBRATE_POINTING#ON_SOURCE. Ignoring this SPW ...')
                            #else:
                            refMaxbw = thebw

                print('Min Freq in DIFFGAIN SPWs: ', minFreq)
                print('Max Freq in DIFFGAIN SPWs: ', maxFreq)
                print('   Max/Min ratio: ', maxFreq/minFreq)
                print('Max BW in DIFFGAIN SIGNAL SPWs:    ', sigMaxbw)
                print('Max BW in DIFFGAIN REFERENCE SPWs: ', refMaxbw)
                print('   Ref/Sig Max BW ratio: ', refMaxbw/sigMaxbw)

                if maxFreq > 1.6611* minFreq: # 1.661 is the largest possible ratio within one band
                    print('The SPWs for intent CALIBRATE_DIFFGAIN fall into different bands. This is an observation with band-to-band phase transfer.')
                    casalog.post('This is an observation with band-to-band phase transfer.', 'INFO')
                    isB2B = True

                elif refMaxbw > 1.5*sigMaxbw:
                    print('The SPWs for CALIBRATE_DIFFGAIN#REFERENCE are wider than for '+diffGainSignalIntent)
                    print('This is an observation with bandwidth switching.')
                    casalog.post('This is an observation with bandwidth switching.', 'INFO')
                    isBWSW = True

                else:
                    print('This observation does not seem to use B2B phase transfer nor BW switching.')

            else:
                casalog.post('Phasecal scans use different SPWs than the scans with intent CALIBRATE_DIFFGAIN#REFERENCE!', 'WARN')
                casalog.post('Will try to treat this as an observation with bandwidth switching but it may fail ...', 'WARN')
                isBWSW = True

        elif (spwsDiffgainsig != []):
            casalog.post('Found intents '+diffGainSignalIntent+' but not CALIBRATE_DIFFGAIN#REFERENCE.', 'WARN')
        elif (spwsDiffgainref != []):
            casalog.post('Found intents CALIBRATE_DIFFGAIN#REFERENCE but not CALIBRATE_DIFFGAIN#SIGNAL nor CALIBRATE_DIFFGAIN#ON_SOURCE.', 'WARN')
    else:
        casalog.post('Found no SPWs for intent CALIBRATE_PHASE#ON_SOURCE !', 'WARN')



    return isB2B, isBWSW
            

#######################################

def listOfIntentsWithFields(msName):
    """
    Return list of intents with field names and ids as a string.
    """
    rval = ''

    intentSources = sfsdr.getIntentsAndSourceNames(msName)

    for k in sorted(intentSources.keys()):
        rval += '# '+k+': '+', '.join(dict.fromkeys(intentSources[k]['name']).keys()) \
                +'\n#     field ids: '+str(intentSources[k]['id'])+'\n'
    
    return rval
    
#######################################


def doQa2ReportGeneration(msName, refant='', isB2B=False, isBWSW=False, iHaveSplitMyScienceSpw=False):
    """
    Generate code which calls the QA2 report generator
    aU.stuffForScienceDataReduction.generateQA2Report()

    """

    print('\n*** doQa2ReportGeneration ***')

    if iHaveSplitMyScienceSpw:
        casalog.post('Reindexing (due to split) of the SPW IDs for the generation of the generateQA2Report call not yet implemented.','WARN')

    print('Gathering information ...')

    casaCmd = "if applyonly != True: es.generateQA2Report(\'"+msName+"\'"

    pardict = {}

    if refant != '':
        pardict['refAnt'] = str(refant)

    if isB2B:

        intentSources = sfsdr.getIntentsAndSourceNames(msName)
        myspwinfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET', caching=True)

        pardict['ant_amp_temporal_caltable'] = msName+'.split.ampli_inf'

        pardict['ant_phase_temporal_caltable'] = msName+'.split.phasehigh_int'

        pardict['ampcalspw'] = ','.join(str(i) for i in myspwinfo.keys())  # the HF SPWs

        pardict['wvrspw'] = str(list(myspwinfo.keys())[0])  # the first of the HF SPWs
        
        pardict['phase_cal'] = intentSources['CALIBRATE_DIFFGAIN']['idstring'][0] # actually the DGC

        pardict['phasecalspw'] = pardict['ampcalspw']

        pardict['checksourceAnalysis'] = False


    elif isBWSW:

        intentSources = sfsdr.getIntentsAndSourceNames(msName)
        myspwinfo = sfsdr.getSpwInfo(msName, intent='OBSERVE_TARGET', caching=True)
        myspwinfoamp = sfsdr.getSpwInfo(msName, intent='CALIBRATE_DIFFGAIN#REFERENCE', caching=True)

        pardict['phase_cal'] = intentSources['CALIBRATE_PHASE']['idstring'][0]

        pardict['target'] = intentSources['OBSERVE_TARGET']['idstring'][0]
 
        pardict['dospw'] = ','.join(str(i) for i in myspwinfo.keys())

        pardict['ampcalspw'] =  ','.join(str(i) for i in myspwinfoamp.keys())

        pardict['wvrspw'] = str(list(myspwinfoamp.keys())[0])  # the first of the SPWs
        
        pardict['phasecalspw'] = pardict['ampcalspw']

    print('Writing code ...')

    if pardict != {}:

        for mypar in pardict.keys():
            if type(mypar) != str:
                casalog.post('Internal ERROR: alls keys of pardict must be strings.', 'SEVERE')
                return casaCmd

            casaCmd += ",\n"

            myparval = pardict[mypar]
            if type(myparval) == str:
                myparval = "'"+myparval+"'"
            else:
                myparval = str(myparval)

            casaCmd += "                                           "+mypar+" = "+myparval

    casaCmd += ")\n"

    return casaCmd

