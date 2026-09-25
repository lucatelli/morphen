# This file contains functions that have been replaced in
# analysisUtils by newer functions which use the ASDM bindings
# library.  We need to keep these for usage on machines that do not
# have this library installed.  - Todd Hunter
from __future__ import print_function  # prevents adding old-style print statements
import os, sys
import math
import numpy as np
from xml.dom import minidom
mylocals = locals()  # used by createCasaTool

casaVersion = None
try:
    import casalith
    casaVersion = casalith.version_string()
except:
    try:
        import casashell # modular casa
        casaVersion = casashell.version_string()
    except:
      # either we are importing into python, or CASA < 6
      if (os.getenv('CASAPATH') is not None):
        import casadef
        if casadef.casa_version >= '5.0.0':
            import casa as mycasa
            if 'cutool' in dir(mycasa):
                cu = mycasa.cutool()
                casaVersion = '.'.join([str(i) for i in cu.version()[:-1]]) + '-' + str(cu.version()[-1])
            else:
                casaVersion = mycasa.casa['build']['version'].split()[0]
        else:
            casaVersion = casadef.casa_version
        print("casaVersion = ", casaVersion)
      else:
        casaVersion = None
        print("os.getenv('CASAPATH') = ", os.getenv('CASAPATH'))
        print("au_noASDMLibrary: You appear to be importing analysisUtils into python (not CASA). version = ", '.'.join([str(i) for i in sys.version_info[:3]]))
if casaVersion is not None:    
    try:
        from taskinit import *
        print("imported casatasks and tools using taskinit *")
    except:
        if casaVersion >= '5.9.9':
            from casatools import quanta as qatool
SECONDS_PER_YEAR = 31556925.

def readCalPointingTable_minidom(sdmfile):
    xmlscans = minidom.parse(sdmfile+'/CalPointing.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    myqa = qatool()
    for rownode in rowlist:
        scandict[fid] = {}
        rowAntennaName = rownode.getElementsByTagName("antennaName")
        antenna = str(rowAntennaName[0].childNodes[0].nodeValue).strip()

        rowDirection = rownode.getElementsByTagName("direction")
#        for r in range(len(rowDirection)):
#            for c in range(len(rowDirection[r].childNodes)):
#                print "%d %d = " % (r,c), rowDirection[r].childNodes[c].nodeValue
        tokens = rowDirection[0].childNodes[0].nodeValue.split()
        azimuth = float(tokens[2])
        elevation = float(tokens[3])

        rowFrequency = rownode.getElementsByTagName("frequencyRange")
        tokens = rowFrequency[0].childNodes[0].nodeValue.split()
        frequency1 = float(tokens[2])
        frequency2 = float(tokens[3])
        frequency = 0.5*(frequency1+frequency2)

        rowRelative = rownode.getElementsByTagName("collOffsetRelative")
        tokens = rowRelative[0].childNodes[0].nodeValue.split()
        azOffset = float(tokens[3])
        elOffset = float(tokens[4])
        azOffset2 = float(tokens[5])
        elOffset2 = float(tokens[6])
        
        rowRelative = rownode.getElementsByTagName("collError")
        tokens = rowRelative[0].childNodes[0].nodeValue.split()
        azError = float(tokens[3])
        elError = float(tokens[4])
        azError2 = float(tokens[5])
        elError2 = float(tokens[6])
        
        rowCalDataId = rownode.getElementsByTagName("calDataId")
        calDataId = str(rowCalDataId[0].childNodes[0].nodeValue)
        scan = int(calDataId.split('_')[1])

        rowpol = rownode.getElementsByTagName("polarizationTypes")
        tokens = rowpol[0].childNodes[0].nodeValue.split()
        poltypes = []
        poltypes.append(str(tokens[2]))
        poltypes.append(str(tokens[3]))

        # start and end times in mjd ns
        rowstart = rownode.getElementsByTagName("startValidTime")
        start = int(rowstart[0].childNodes[0].nodeValue)/1000000000
        startmjd = start/86400.0
        t = myqa.quantity(startmjd,'d')
        starttime = call_qa_time(t,form="ymd",prec=8)
        rowend = rownode.getElementsByTagName("endValidTime")
        end = int(rowend[0].childNodes[0].nodeValue)
        endmjd = float(end)*1.0E-9/86400.0
        t = myqa.quantity(endmjd,'d')
        endtime = call_qa_time(t,form="ymd",prec=8)

        scandict[fid]['startValidTime'] = start
        scandict[fid]['endValidTime'] = end
        scandict[fid]['start'] = starttime
        scandict[fid]['end'] = endtime
        scandict[fid]['startmjd'] = startmjd
        scandict[fid]['endmjd'] = endmjd
        scandict[fid]['startmjdsec'] = startmjd*86400
        scandict[fid]['endmjdsec'] = endmjd*86400
        timestr = starttime+'~'+endtime
        scandict[fid]['azimuth'] = azimuth
        scandict[fid]['elevation'] = elevation
        scandict[fid]['antenna'] = antenna
        scandict[fid]['frequency'] = frequency
        scandict[fid]['azOffset'] = azOffset
        scandict[fid]['elOffset'] = elOffset
        scandict[fid]['azOffset2'] = azOffset2
        scandict[fid]['elOffset2'] = elOffset2
        scandict[fid]['azError'] = azError
        scandict[fid]['elError'] = elError
        scandict[fid]['azError2'] = azError2
        scandict[fid]['elError2'] = elError2
        scandict[fid]['scan'] = scan
        scandict[fid]['duration'] = (endmjd-startmjd)*86400
        scandict[fid]['polarizationTypes'] = poltypes
        fid += 1
    print('  Found ',rowlist.length,' rows in CalPointing.xml')
    myqa.done()
    # return the dictionary for later use
    return scandict

def readCalDataFromASDM_minidom(asdm):
    """
    Returns a dictionary, keyed by 'scan' and 'calDataId' where the value of scan is a tuple
    """
    xmlscans = minidom.parse(asdm+'/CalData.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    for fid,rownode in enumerate(rowlist):
        scandict[fid] = {}
        rowscan = rownode.getElementsByTagName("scanSet")
        tokens = rowscan[0].childNodes[0].nodeValue.split()
        nscans = int(tokens[1])
        scan = []
        for i in range(nscans):
            scan.append(int(tokens[2+i]))
        rowcaldataid = rownode.getElementsByTagName("calDataId")
        caldataid = str(rowcaldataid[0].childNodes[0].nodeValue).strip()
        scandict[fid]['calDataId'] = caldataid
        scandict[fid]['scan'] = tuple(scan)
    return(scandict)

def getObservatoryNameFromASDM_minidom(asdm):
    execblock = asdm + '/ExecBlock.xml'
    if not os.path.exists(execblock):
        print("Could not open %s" % (execblock))
        return
    xmlscans = minidom.parse(execblock)
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    row = rowlist[0].getElementsByTagName("telescopeName")
    myName = str(row[0].childNodes[0].nodeValue).strip(' ')
    return myName

def getRADecForFieldFromASDM_minidom(asdm, field, ephemerisFields):
    xmlfields = minidom.parse(asdm+'/Field.xml')
    rowlist = xmlfields.getElementsByTagName("row")
    fields = []
    ra = None; dec = None
    for rownode in rowlist:
        rownumLO = rownode.getElementsByTagName("fieldId")
        fieldId = int(str(rownumLO[0].childNodes[0].nodeValue).split('_')[1])
        fields.append(fieldId)
        if (fieldId == field):
            if field in ephemerisFields:
                radec = ephemerisFields[field]['direction']
                ra = radec # leave dec as None to signify sexagesimal string
            else:
                rownumLO = rownode.getElementsByTagName("phaseDir")
                ra = float(str(rownumLO[0].childNodes[0].nodeValue).strip().split(' ')[3])
                dec = float(str(rownumLO[0].childNodes[0].nodeValue).strip().split(' ')[4])
            break
    return ra, dec

def qa0flags_minidom(asdm):
    flag = asdm + '/Flag.xml'
    xmlscans = minidom.parse(flag)
    rowlist = xmlscans.getElementsByTagName("row")
    reasons = []; commands = []; durations = []
    tsys = 0; trx = 0
    tsysAntennas = []; tsysBasebands = []; tsysPols = []
    trxAntennas = []; trxBasebands = []; trxPols = []
    for i,rownode in enumerate(rowlist):
        rowStart = rownode.getElementsByTagName("startTime")
        rowEnd = rownode.getElementsByTagName("endTime")
        rowReason = rownode.getElementsByTagName("reason")
        reasons.append(str(rowReason[0].childNodes[0].nodeValue).strip().replace(' ','_'))
        if reasons[i].find('AFD05') > 0:
            if reasons[i].find('TSYS') > 0:
                tsys += 1
                tsysAntennas.append(reasons[i].split('Ant_')[1].split('_')[0])
                tsysBasebands.append(int(reasons[i].split('(BB_')[1].split(')')[0]))
                tsysPols.append(reasons[i].split('_Pol_')[1].split('_')[0])
            elif reasons[i].find('TRX') > 0:
                trx += 1
                trxAntennas.append(reasons[i].split('Ant_')[1].split('_')[0])
                trxBasebands.append(int(reasons[i].split('(BB_')[1].split(')')[0]))
                trxPols.append(reasons[i].split('_Pol_')[1].split('_')[0])
        startTime = int(rowStart[0].childNodes[0].nodeValue)
        endTime =   int(rowEnd[0].childNodes[0].nodeValue)
        commands.append([startTime,endTime])
    return reasons, commands, tsys, trx, tsysAntennas, tsysBasebands, tsysPols, trxAntennas, trxBasebands, trxPols

def getTsysFromSysCal_minidom(asdm, sqld, channel, scanlist, verbose):
    xmlscans = minidom.parse(asdm+'/SpectralWindow.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    wvrSpws = []
    baseband = {}
    for rownode in rowlist:
        name = rownode.getElementsByTagName("name")
        name = str(name[0].childNodes[0].nodeValue)
        spwID = rownode.getElementsByTagName("spectralWindowId")
        spw = int(str(spwID[0].childNodes[0].nodeValue).split('_')[1])
        basebandName = rownode.getElementsByTagName("basebandName")
        basebandName = str(basebandName[0].childNodes[0].nodeValue).split('_')  # e.g. "BB_1"
        if (len(basebandName) < 2):
            baseband[spw] = -1  # NOBB
        else:
            baseband[spw] = int(basebandName[1])
        if (name.find('WVR') >= 0):
            wvrSpws.append(spw)
    xmlscans = minidom.parse(sqld)
    rowlist = xmlscans.getElementsByTagName("row")
    scandict = {}
    duplicates = 0
    nchans = []
    if len(rowlist) == 0:
        print("No rows found in SysCal.xml, probably because all the data are now in the binary file.")
        return
    for rownode in rowlist:
        antennaID = rownode.getElementsByTagName("antennaId")
        antenna = int(str(antennaID[0].childNodes[0].nodeValue).split('_')[1])
        if (antenna not in list(scandict.keys())):
            scandict[antenna] = {}
        spwID = rownode.getElementsByTagName("spectralWindowId")
        asdmspw = int(str(spwID[0].childNodes[0].nodeValue).split('_')[1])
        subtract = len(np.where(asdmspw > np.array(wvrSpws))[0])-1 
        spw = asdmspw-subtract # translate this to actual spw number, as only 1 WVR spw is real
        if (spw not in list(scandict[antenna].keys())):
            scandict[antenna][spw] = {}
        scandict[antenna][spw]['baseband'] = baseband[asdmspw]
        scandict[antenna][spw]['scans'] = {}
        timeData = rownode.getElementsByTagName("timeInterval")
        timeStamp = int(str(timeData[0].childNodes[0].nodeValue).split()[0])
        timeInterval = int(str(timeData[0].childNodes[0].nodeValue).split()[1])
        timeCenter = timeStamp-timeInterval/2
        timeCenterMJD = timeCenter*1e-9/86400.
        scan = -1
        for s in list(scanlist.keys()):
            if (scanlist[s]['endmjd'] >= timeCenterMJD and scanlist[s]['startmjd'] <= timeCenterMJD):
                scan = s
        if (scan not in list(scandict[antenna][spw]['scans'].keys())):
            scandict[antenna][spw]['scans'][scan] = {}
        tsysSpectrum = rownode.getElementsByTagName("tsysSpectrum")
        values = tsysSpectrum[0].childNodes[0].nodeValue.split()
        npol = int(values[1])  # or is it [0] ?
        nchan = int(values[2])
        nchans.append(nchan)
        for pol in range(npol):
            tsys = []
            for i in range(pol*nchan, nchan*(pol+1)):
                tsys.append(float(values[3+i]))
            if (channel is None):
                scandict[antenna][spw]['scans'][scan][pol] = tsys
            else:
                scandict[antenna][spw]['scans'][scan][pol] = [tsys[channel]]
    if verbose:
        print('  Found ',rowlist.length,' Tsys rows in SysCal.xml (median # chans = %.0f)' % (np.median(nchans)))
    return scandict

def getPolarizationsFromASDM_minidom(asdm):
    xmlscans = minidom.parse(asdm+'/Polarization.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    mydict = {}
    for rownode in rowlist:
        rowscan = rownode.getElementsByTagName("numCorr")
        tokens = rowscan[0].childNodes[0].nodeValue.split()
        numCorr = int(tokens[0])
        rowscan = rownode.getElementsByTagName("corrType")
        tokens = rowscan[0].childNodes[0].nodeValue.split()
        corrTypes = [str(i) for i in tokens[2:]]
        rowcaldataid = rownode.getElementsByTagName("polarizationId")
        polarizationId = int(str(rowcaldataid[0].childNodes[0].nodeValue).split('_')[1])
        mydict[polarizationId] = {'numCorr': numCorr, 'corrTypes': corrTypes}
    return mydict

def getScienceSpwTransitionsFromASDM_minidom(asdm):
    xmlscans = minidom.parse(asdm+'/Source.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    transitions = {}
    scienceSpwsASDM = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("velRefCode")
        if (len(row) > 0):
            velRefCode = str(row[0].childNodes[0].nodeValue)
            row = rownode.getElementsByTagName("spectralWindowId")
            spw = int(str(row[0].childNodes[0].nodeValue).split('_')[1])
            scienceSpwsASDM.append(spw)
            row = rownode.getElementsByTagName("transition")
            transition = str(row[0].childNodes[0].nodeValue).split()[2].strip('"')
            transitions[spw] = transition
    return transitions, scienceSpwsASDM

def readPointingModelFromASDM_minidom(asdm, antennaNames):
    xmlscans = minidom.parse(asdm+'/PointingModel.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    coeff = {}
    for rownode in rowlist:
        rowAntennaId = rownode.getElementsByTagName("antennaId")
        antennaId = str(rowAntennaId[0].childNodes[0].nodeValue)
        rowCoeffName = rownode.getElementsByTagName("coeffName")
        coeffName = str(rowCoeffName[0].childNodes[0].nodeValue)
        coeffNames = ', '.join(coeffName.split()[2:]).replace('"','')
        rowAssocNature = rownode.getElementsByTagName("assocNature")
        assocNature = str(rowAssocNature[0].childNodes[0].nodeValue).strip()
        rowPolarizationType = rownode.getElementsByTagName("polarizationType")
        polarizationType = str(rowPolarizationType[0].childNodes[0].nodeValue)
        rowCoeffVal = rownode.getElementsByTagName("coeffVal")
        coeffVal = str(rowCoeffVal[0].childNodes[0].nodeValue)
        tokens = coeffVal.split()
        antennaId = int(antennaId.split('_')[1])
        antennaName = antennaNames[antennaId]
        if (antennaName not in list(coeff.keys())):
            coeff[antennaName] = {}
        if (polarizationType not in list(coeff[antennaName].keys())):
            coeff[antennaName][polarizationType] = {}
        if (assocNature not in list(coeff[antennaName][polarizationType].keys())):
            coeff[antennaName][polarizationType][assocNature] = []
        coeff[antennaName][polarizationType][assocNature].append([float(token) for token in tokens[2:]])
        coeff[antennaName]['id'] = antennaId
    return coeff,coeffNames

def readSysCal_minidom(asdm):
    xmlscans = minidom.parse(asdm+'/SysCal.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    duplicates = 0
    scandict = {}
    for rownode in rowlist:
        antennaID = rownode.getElementsByTagName("antennaId")
        antenna = int(str(antennaID[0].childNodes[0].nodeValue).split('_')[1])
        if (antenna not in list(scandict.keys())):
            scandict[antenna] = {}
        spwID = rownode.getElementsByTagName("spectralWindowId")
        spw = int(str(spwID[0].childNodes[0].nodeValue).split('_')[1])
        if (spw not in list(scandict[antenna].keys())):
            scandict[antenna][spw] = []
        timeData = rownode.getElementsByTagName("timeInterval")
        timeStamp = int(str(timeData[0].childNodes[0].nodeValue).split()[0])
        timeInterval = int(str(timeData[0].childNodes[0].nodeValue).split()[1])
        if (timeStamp-timeInterval/2 in scandict[antenna][spw]):
            print("Duplicate seen!")
            duplicates += 1
        else:
            scandict[antenna][spw].append(timeStamp-timeInterval/2)
        fid += 1
    print('  Found ',rowlist.length,' Tsys rows in SysCal.xml')
    if rowlist.length > 0:
        print("%d duplicates found" % (duplicates))
    else:
        print("But this is likely because the information is stored in SysCal.bin rather than SysCal.xml, in which case you need to use the ASDM python libraries.")
    return scandict

def getIntegrationTimeFromASDM_minidom(asdm):
    xmlscans = minidom.parse(asdm+'/Main.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    integrationTimes = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("numIntegration")
        tokens = row[0].childNodes[0].nodeValue.split()
        numIntegration = int(str(tokens[0]))
        row = rownode.getElementsByTagName("interval")
        tokens = row[0].childNodes[0].nodeValue.split()
        interval = float(str(tokens[0])) * 1e-9
        integrationTimes.append(interval / numIntegration)
    return integrationTimes

def readFeedFromASDM_minidom(asdm):
    xmlscans = minidom.parse(asdm+'/Feed.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    values0 = []
    values1 = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("timeInterval")
        tokens = row[0].childNodes[0].nodeValue.split()
        values0.append(float(tokens[0])*1e-9)
        values1.append(float(tokens[1])*1e-9)
    values0 = np.unique(values0)
    values1 = np.unique(values1)
    return values0, values1

def representativeFrequencyFromASDM_minidom(asdm):
    freq = None
    f = open(asdm+'/SBSummary.xml')
    for line in f.readlines():
        loc = line.find('representativeFrequency')
        if (loc > 0):
            tokens = line[loc:].split()
            freqArg = tokens[2]+tokens[3]
        loc = line.find('representativeWindow')
        if (loc > 0):
            tokens = line[loc:].split()
            spw = str(tokens[2]).strip()
    f.close()
    return freqArg, spw

def getLOsFromASDM_minidom(sdmfile):
    xmlscans = minidom.parse(sdmfile+'/Receiver.xml')
    fid = 0
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    for rownode in rowlist:
        scandict[fid] = {}
        rownumLO = rownode.getElementsByTagName("numLO")
        numLO = int(rownumLO[0].childNodes[0].nodeValue)
        rowfreqLO = rownode.getElementsByTagName("freqLO")
        rowreceiverId = rownode.getElementsByTagName("receiverId")
        receiverId = int(rowreceiverId[0].childNodes[0].nodeValue)
        rowfrequencyBand = rownode.getElementsByTagName("frequencyBand")
        frequencyBand = str(rowfrequencyBand[0].childNodes[0].nodeValue)
        freqLO = []
        r = list(filter(None,(rowfreqLO[0].childNodes[0].nodeValue).split(' ')))
        for i in range(2,len(r)):
            freqLO.append(float(r[i]))
        
        rowspwid = rownode.getElementsByTagName("spectralWindowId")
        spwid = int(str(rowspwid[0].childNodes[0].nodeValue).split('_')[1])
        scandict[fid]['spectralWindowId'] = spwid
        scandict[fid]['freqLO'] = freqLO
        scandict[fid]['numLO'] = numLO
        scandict[fid]['receiverId'] = receiverId
        scandict[fid]['frequencyBand'] = frequencyBand
        fid +=1
    return scandict

def requestedResolutionFromASDM_minidom(asdm):
    f = open(asdm+'/SBSummary.xml')
    minAcceptableResolution = 0
    maxAcceptableResolution = 0
    for line in f.readlines():
        loc = line.find('minAcceptableAngResolution')
        if (loc >= 0):
            tokens = line[loc:].split()
            minAcceptableResolution = float(tokens[2])
            minUnits = tokens[3]
        loc = line.find('maxAcceptableAngResolution')
        if (loc >= 0):
            tokens = line[loc:].split()
            maxAcceptableResolution = float(tokens[2])
            maxUnits = tokens[3]
    f.close()
    return minAcceptableResolution, maxAcceptableResolution

def restFrequenciesASDM_minidom(asdm):
    restFreqs = {}
    xmlscans = minidom.parse(asdm+'/Source.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    for rownode in rowlist:
        row = rownode.getElementsByTagName("velRefCode")
        if (len(row) > 0):
            row = rownode.getElementsByTagName("spectralWindowId")
            spw = int(str(row[0].childNodes[0].nodeValue).split('_')[1])
            row = rownode.getElementsByTagName("restFrequency")
            tokens = (row[0].childNodes[0].nodeValue).split()
            restFreqs[spw] = []
            for i in range(int(tokens[1])):
                restFreqs[spw].append(float(tokens[2+i]))
    return restFreqs

def properMotionASDM_minidom(asdm, arcsecPerYear):
    xmlscans = minidom.parse(asdm+'/Source.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    properMotion = {}
    for rownode in rowlist:
        row = rownode.getElementsByTagName("sourceName")
        tokens = (row[0].childNodes[0].nodeValue).split()
        sourceName = str(tokens[0])
        row = rownode.getElementsByTagName("properMotion")
        if (len(row) > 0 and sourceName not in properMotion):
            tokens = (row[0].childNodes[0].nodeValue).split()
            properMotion[sourceName] = [float(tokens[2]), float(tokens[3])]
            if arcsecPerYear:
                properMotion[sourceName] = list(np.degrees(np.array(properMotion[sourceName]))*3600*SECONDS_PER_YEAR)
    return properMotion

def readWeatherFromASDM_minidom(sdmfile):
    xmlscans = minidom.parse(sdmfile+'/Weather.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    for rownode in rowlist:
        scandict[fid] = {}
        row = rownode.getElementsByTagName("timeInterval")
        tokens = row[0].childNodes[0].nodeValue.split()
        scandict[fid]['timeInterval'] = float(tokens[0])*1e-9  # MJD seconds
        row = rownode.getElementsByTagName("pressure")
        scandict[fid]['pressure'] = float(row[0].childNodes[0].nodeValue)*0.01 # mbar
        row = rownode.getElementsByTagName("relHumidity")
        scandict[fid]['relHumidity'] = float(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("temperature")
        scandict[fid]['temperature'] = float(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("windDirection")
        scandict[fid]['windDirection'] = float(row[0].childNodes[0].nodeValue)*180/math.pi  # degrees
        row = rownode.getElementsByTagName("windSpeed")
        scandict[fid]['windSpeed'] = float(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("windMax")
        scandict[fid]['windMax'] = float(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("stationId")
        scandict[fid]['stationId'] = int((str(row[0].childNodes[0].nodeValue)).split('_')[1])
        fid += 1
    return scandict

def getCorrelatorName_minidom(asdm):
    execblock = asdm + '/CorrelatorMode.xml'
    if (os.path.exists(execblock) == False):
        print("Could not open %s" % (execblock))
        return
    xmlscans = minidom.parse(execblock)
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    row = rowlist[0].getElementsByTagName("correlatorName")
    myName = str(row[0].childNodes[0].nodeValue)
    return myName

def getEphemerisFieldsFromASDM_minidom(asdm, keyBy, verbose):        
    xmlscans = minidom.parse(asdm+'/Ephemeris.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    ephemeris = {}
    for rownode in rowlist:
        myrow = rownode.getElementsByTagName("ephemerisId")
        ephemerisId = int(myrow[0].childNodes[0].nodeValue)
        if ephemerisId not in ephemeris:
            myrow = rownode.getElementsByTagName("dir")
            rad = [float(myrow[0].childNodes[0].nodeValue.split()[3]), float(myrow[0].childNodes[0].nodeValue.split()[4])]
#            radec = rad2radec(rad, verbose=False)
            myrow = rownode.getElementsByTagName("timeInterval")
            mjd = float(myrow[0].childNodes[0].nodeValue.split()[0])*1e-9/86400.
            ephemeris[ephemerisId] = {'direction': rad, 'mjd': mjd}
    xmlscans = minidom.parse(asdm+'/Field.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    ephemerisByName = {}
    ephemerisByFieldID = {}
    for rownode in rowlist:
        myrow = rownode.getElementsByTagName("ephemerisId")
        if len(myrow) > 0:
            ephemerisId = int(myrow[0].childNodes[0].nodeValue)
            myrow = rownode.getElementsByTagName("fieldName")
            fieldName = str(myrow[0].childNodes[0].nodeValue).strip()
            ephemerisByName[fieldName] = ephemeris[ephemerisId]
            myrow = rownode.getElementsByTagName("fieldId")
            fieldId = int(str(myrow[0].childNodes[0].nodeValue).split('_')[1])
            ephemerisByFieldID[fieldId] = ephemeris[ephemerisId]
    return ephemeris, ephemerisByFieldID, ephemerisByName

def getSpwsFromASDM_minidom(sdmfile, minnumchan=0, dropExtraWVRSpws=True):
    fid = 0
    firstWVR = True
    scandict = {}
    xmlscans = minidom.parse(sdmfile+'/SpectralWindow.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    for rowNumber, rownode in enumerate(rowlist):
        rowspwid = rownode.getElementsByTagName("spectralWindowId")
        spwid = int(str(rowspwid[0].childNodes[0].nodeValue).split('_')[1])
        rownumLO = rownode.getElementsByTagName("numChan")
        numChan = int(rownumLO[0].childNodes[0].nodeValue)
        if (numChan < minnumchan): 
            if (numChan == 4):
                if (firstWVR):
                    fid += 1
                    firstWVR = False
            else:
                fid += 1
            continue
        scandict[fid] = {}
        rownumLO = rownode.getElementsByTagName("refFreq")
        refFreq = float(rownumLO[0].childNodes[0].nodeValue)
        rownumNOBB = rownode.getElementsByTagName("basebandName")
        noBB = rownumNOBB[0].childNodes[0].nodeValue
        rownumName = rownode.getElementsByTagName("name")
        name = rownumName[0].childNodes[0].nodeValue.strip()
        rownumWF = rownode.getElementsByTagName("windowFunction")
        windowFunction = str(rownumWF[0].childNodes[0].nodeValue)
        try:
            rownumWF = rownode.getElementsByTagName("effectiveBw")
            effectiveBw = float(rownumWF[0].childNodes[0].nodeValue)
        except:
            rownumWF = rownode.getElementsByTagName("effectiveBwArray")
            effectiveBw = float(list(filter(None,(rownumWF[0].childNodes[0].nodeValue).split()))[2])
        try:
            rownumWF = rownode.getElementsByTagName("chanWidth")
            chanWidth = float(rownumWF[0].childNodes[0].nodeValue)
        except:
            rownumWF = rownode.getElementsByTagName("chanWidthArray")
            chanWidth = float(list(filter(None,(rownumWF[0].childNodes[0].nodeValue).split()))[2])
        try:
            rownumWF = rownode.getElementsByTagName("resolution")
            resolution = float(rownumWF[0].childNodes[0].nodeValue)
        except:
            rownumWF = rownode.getElementsByTagName("resolutionArray")
            resolution = float(list(filter(None,(rownumWF[0].childNodes[0].nodeValue).split()))[2])
        scandict[fid]['spectralWindowId'] = spwid
        scandict[fid]['numChan'] = numChan
        scandict[fid]['refFreq'] = refFreq
        scandict[fid]['windowFunction'] = windowFunction
        scandict[fid]['effectiveBw'] = effectiveBw
        scandict[fid]['resolution'] = resolution
        scandict[fid]['chanWidth'] = chanWidth
        scandict[fid]['name'] = name
        if (noBB == 'NOBB'):
            scandict[fid]['basebandNumber'] = 0
        else:
            scandict[fid]['basebandNumber'] = int(noBB.split('_')[-1])
        try:
            rownumLO = rownode.getElementsByTagName("chanFreqStart")
            chanFreqStart = float(rownumLO[0].childNodes[0].nodeValue)
            rownumLO = rownode.getElementsByTagName("chanFreqStep")
            chanFreqStep = float(rownumLO[0].childNodes[0].nodeValue)
            if ((numChan % 2) == 1):
                centerFreq = chanFreqStart + (numChan/2)*chanFreqStep
            else:
                centerFreq = chanFreqStart + (numChan-1)*0.5*chanFreqStep
        except:
            try:
                rownumLO = rownode.getElementsByTagName("chanFreqArray")
                r = list(filter(None,(rownumLO[0].childNodes[0].nodeValue).split(' ')))
                freqLO = []
                for i in range(2,len(r)):
                    freqLO.append(float(r[i]))
                centerFreq = np.mean(freqLO)
                chanFreqStart = freqLO[0]
            except:
                print("Did not find chanFreqStart nor chanFreqArray on row=%d, spw=%d" % (fid,spwid))
                scandict[fid]['centerFreq'] = 0
                continue
        scandict[fid]['centerFreq'] = centerFreq
        scandict[fid]['chanFreqStart'] = chanFreqStart
        if (refFreq > centerFreq):
            scandict[fid]['sideband'] = -1
        else:
            scandict[fid]['sideband'] = +1
        if (scandict[fid]['numChan'] != 4 or firstWVR or dropExtraWVRSpws==False):
            if (firstWVR and scandict[fid]['numChan'] == 4):
                scandict[fid]['sideband'] = 0
                firstWVR = False
            fid += 1
    return scandict
    
def getScienceSpwsFromASDM_minidom(asdm, returnTransitions=False):
    xmlscans = minidom.parse(asdm+'/Source.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    scienceSpwsASDM = []
    transitions = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("velRefCode")
        if (len(row) > 0):
            velRefCode = str(row[0].childNodes[0].nodeValue)
            row = rownode.getElementsByTagName("spectralWindowId")
            spw = int(str(row[0].childNodes[0].nodeValue).split('_')[1])
            scienceSpwsASDM.append(spw)
            row = rownode.getElementsByTagName("transition")
            transitions.append(str(row[0].childNodes[0].nodeValue))
    if returnTransitions:
        return scienceSpwsASDM, transitions
    else:
        return scienceSpwsASDM

def getFieldsFromASDM_minidom(asdm):
    mydict = {}
    mydict2 = {}
    xmlscans = minidom.parse(asdm+'/Field.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    for rownode in rowlist:
        rowpwv = rownode.getElementsByTagName("fieldId")
        fieldid = int(rowpwv[0].childNodes[0].nodeValue.split('_')[1])
        rowpwv = rownode.getElementsByTagName("fieldName")
        fieldname = str(rowpwv[0].childNodes[0].nodeValue).strip()
        mydict[fieldid] = fieldname
        mydict2[fieldname] = fieldid
    return mydict, mydict2

def getObservationEndDateFromASDM_minidom(asdm):
    """
    Returns the end date/time and MJD seconds in the specified ASDM.
    --  based on getObservationStartDateFromASDM
    """
#    uid___A002_X54d35d_X761/ExecBlock.xml:    <startTime>ALMA</startTime>
    execblock = asdm + '/ExecBlock.xml'
    if (os.path.exists(execblock) == False):
        print("Could not open %s" % (execblock))
        return
    xmlscans = minidom.parse(execblock)
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    row = rowlist[0].getElementsByTagName("endTime")
    mjdsec = int(row[0].childNodes[0].nodeValue) * 1e-9
    return(mjdsec)

def call_qa_time(arg, form='', prec=0, showform=False):
    """
    This is a wrapper for qa.time(), which in casa 4.0.0 returns a list 
    of strings instead of just a scalar string.  
    arg: a time quantity
    - Todd Hunter
    """
    if (type(arg) == dict):
        if (type(arg['value']) == list or 
            type(arg['value']) == np.ndarray):
            if (len(arg['value']) > 1):
                print("WARNING: call_qa_time() received a dictionary containing a list of length=%d rather than a scalar. Using first value." % (len(arg['value'])))
            arg['value'] = arg['value'][0]
    myqa = qatool()
    result = myqa.time(arg, form=form, prec=prec, showform=showform)
    myqa.done()
    if (type(result) == list or type(result) == np.ndarray):
        return(result[0])
    else:
        return(result)

def readwvr_minidom(sdmfile, verbose=False):
    """
    This function reads the CalWVR.xml table from the ASDM and returns a
    dictionary containing: 'start', 'end', 'startmjd', 'endmjd',
    'startmjdsec', 'endmjdsec',
    'timerange', 'antenna', 'water', 'duration'.
    'water' is the zenith PWV in meters.
    This function is called by au.readwvr(). -- Todd Hunter
    """
    if (not os.path.exists(sdmfile)):
        print("readwvr_minidom(): Could not find file = ", sdmfile)
        return
    xmlscans = minidom.parse(sdmfile+'/CalWVR.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    myqa = createCasaTool(qatool) # qatool()
    pathDict = {}
    for rownode in rowlist:
        rowpwv = rownode.getElementsByTagName("water")
        pwv = float(rowpwv[0].childNodes[0].nodeValue)
        water = pwv
        scandict[fid] = {}

        # start and end times in mjd ns
        rowstart = rownode.getElementsByTagName("startValidTime")
        start = int(rowstart[0].childNodes[0].nodeValue)
        startmjd = float(start)*1.0E-9/86400.0
        t = myqa.quantity(startmjd,'d')
        starttime = call_qa_time(t,form="ymd",prec=8)
        rowend = rownode.getElementsByTagName("endValidTime")
        end = int(rowend[0].childNodes[0].nodeValue)
        endmjd = float(end)*1.0E-9/86400.0
        t = myqa.quantity(endmjd,'d')
        endtime = call_qa_time(t,form="ymd",prec=8)
        # antenna
        rowantenna = rownode.getElementsByTagName("antennaName")
        antenna = str(rowantenna[0].childNodes[0].nodeValue).strip()
        if antenna not in list(pathDict.keys()):
            pathDict[antenna] = [[],[],[],[]]

        rowfreq = rownode.getElementsByTagName("chanFreq")
        freqs = []
        freqs.append(float(list(filter(None,(rowfreq[0].childNodes[0].nodeValue).split()))[2]))
        freqs.append(float(list(filter(None,(rowfreq[0].childNodes[0].nodeValue).split()))[3]))
        freqs.append(float(list(filter(None,(rowfreq[0].childNodes[0].nodeValue).split()))[4]))
        freqs.append(float(list(filter(None,(rowfreq[0].childNodes[0].nodeValue).split()))[5]))

        rowpath = rownode.getElementsByTagName("pathCoeff")
        pathCoeffs = 4*[0]
        pathCoeffs[0] = float(list(filter(None,(rowpath[0].childNodes[0].nodeValue).split()))[4])
        pathCoeffs[1] = float(list(filter(None,(rowpath[0].childNodes[0].nodeValue).split()))[5])
        pathCoeffs[2] = float(list(filter(None,(rowpath[0].childNodes[0].nodeValue).split()))[6])
        pathCoeffs[3] = float(list(filter(None,(rowpath[0].childNodes[0].nodeValue).split()))[7])
        for i in range(4):
            pathDict[antenna][i].append(pathCoeffs[i])
            
        scandict[fid]['start'] = starttime
        scandict[fid]['end'] = endtime
        scandict[fid]['startmjd'] = startmjd
        scandict[fid]['endmjd'] = endmjd
        scandict[fid]['startmjdsec'] = startmjd*86400
        scandict[fid]['endmjdsec'] = endmjd*86400
        timestr = starttime+'~'+endtime
        scandict[fid]['timerange'] = timestr
        scandict[fid]['antenna'] = antenna
        scandict[fid]['water'] = water
        scandict[fid]['duration'] = (endmjd-startmjd)*86400
        fid += 1

    if verbose: print('  Found ',rowlist.length,' rows in CalWVR.xml')
    myqa.done()
    return scandict, pathDict

def getObservationStartDateFromASDM_minidom(asdm):
    """
    Returns the start date/time and MJD seconds in the specified ASDM.
    example:  ('2016-05-12 07:14:19 UT', 4969754058.985001)
    -- Todd Hunter
    """
#    uid___A002_X54d35d_X761/ExecBlock.xml:    <startTime>ALMA</startTime>
    execblock = asdm + '/ExecBlock.xml'
    if (os.path.exists(execblock) == False):
        print("Could not open %s" % (execblock))
        return
    xmlscans = minidom.parse(execblock)
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    row = rowlist[0].getElementsByTagName("startTime")
    mjdsec = int(row[0].childNodes[0].nodeValue) * 1e-9
    return(mjdsec)

def getScienceSpwsFromASDM_minidom(asdm, returnTransitions=False):
    """
    Gets the spw IDs for the science target spws in an ASDM.  Single-channel spws are
    ignored.  The IDs are those that will appear in a measurement set once it is imported.
    It determines science spws by the presence of the "velRefCode" tag, which does not
    work for the image sidebands of 90deg Walsh SpectralSpecs in Cycle 5 and 6.
    -Todd Hunter
    """
    xmlscans = minidom.parse(asdm+'/Source.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    scienceSpwsASDM = []
    transitions = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("velRefCode")
        if (len(row) > 0):
            velRefCode = str(row[0].childNodes[0].nodeValue)
            row = rownode.getElementsByTagName("spectralWindowId")
            spw = int(str(row[0].childNodes[0].nodeValue).split('_')[1])
            scienceSpwsASDM.append(spw)
            row = rownode.getElementsByTagName("transition")
            transitions.append(str(row[0].childNodes[0].nodeValue))
    scienceSpwsASDM = np.unique(scienceSpwsASDM)
#    print("A) scienceSpws in ASDM: " , scienceSpwsASDM)
    spwmap = asdmspwmap_minidom(asdm)
    scienceSpws = []
    scienceTransitions = []
    mydict = getSpwsFromASDM_minidom(asdm)
    for i,spw in enumerate(scienceSpwsASDM):
        myspw = spwmap.index(spw)
        if mydict[myspw]['numChan'] > 1:
            scienceSpws.append(myspw)
            scienceTransitions.append(transitions[i])
#    print("A) scienceSpws in ms: " , scienceSpws)
    if returnTransitions:
        return scienceSpws, scienceTransitions
    else:
        return scienceSpws

def getIntentsFromASDM_minidom(asdm, stripPrefix=False, byscan=False):
    mydict = {}
    byscandict = {}
    xmlscans = minidom.parse(asdm+'/Scan.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    for rownode in rowlist:
        rowpwv = rownode.getElementsByTagName("fieldName")
        # The following is necessary for the case that the SB was terminated
        # early, in which case the final entry of the Scan.xml might not have a <fieldName>
        if (rowpwv == []): continue
        names = rowpwv[0].childNodes[0].nodeValue.split()
        fieldname = str(names[2]).strip('"')
        if (fieldname not in list(mydict.keys())):
            mydict[fieldname] = []
        rowintent = rownode.getElementsByTagName("scanIntent")
        tokens = rowintent[0].childNodes[0].nodeValue.split()
        numIntents = int(tokens[1])
        scanNumber = int(rownode.getElementsByTagName("scanNumber")[0].childNodes[0].nodeValue.split()[0])
        for i in range(numIntents):
            intent = str(tokens[i+2])
            if stripPrefix:
                intent = intent[intent.find('_')+1:]
            if (intent not in mydict[fieldname]):
                mydict[fieldname].append(intent)
            if (intent not in list(byscandict.keys())):
                byscandict[intent] = []
            byscandict[intent].append(scanNumber)
    return mydict, byscandict

def readAntennasFromASDM_minidom(sdmfile, stations=False, diameters=False, 
                                 verbose=True, useMinidom=False):
    """
    Returns a tuple of 3 lists: antennaNameList, stationList, dishDiameterList
    """
    antList = []
    stationList = []
    dishDiameterList = []
    xmlscans = minidom.parse(sdmfile+'/Antenna.xml')
    rowlist = xmlscans.getElementsByTagName("row")
    for rownode in rowlist:
        row = rownode.getElementsByTagName("name")
        tokens = row[0].childNodes[0].nodeValue.split()
        antList.append(str(tokens[0]))
        if (stations):
            row = rownode.getElementsByTagName("stationId")
            tokens = row[0].childNodes[0].nodeValue.split()
            stationList.append(int(str(tokens[0]).split('_')[1]))
        if (diameters):
            row = rownode.getElementsByTagName("dishDiameter")
            tokens = row[0].childNodes[0].nodeValue.split()
            dishDiameterList.append(float(str(tokens[0])))
    return antList, stationList, dishDiameterList

def readTimesFromAnnotationTable_minidom(asdm):
    print("Observation start date: ", getObservationStartDateFromASDM_minidom(asdm)[0])
    print("Representative frequency: %.3f GHz" % (representativeFrequencyFromASDM_minidom(asdm,verbose=False)))
    xmlscans = minidom.parse(asdm+'/Annotation.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    spws = getScienceSpwsFromASDM_minidom(asdm)
    mydict = {}
    previousIntent = None
    spectralSpec = 0
    for i,rownode in enumerate(rowlist):
        row = rownode.getElementsByTagName("issue")
        issue = str(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("details")
        details = str(row[0].childNodes[0].nodeValue)
        if issue == "Parameter Optimization":
#            print "Parsing line: ", details
            intent, sourcename, originalIntegration, optimumIntegration, adoptedIntegration = details.split(',')
            intent = intent.upper()
            if intent != previousIntent and previousIntent is not None:
                spectralSpec = 0
            previousIntent = intent
            if spectralSpec not in mydict:
                mydict[spectralSpec] = {}
            elif intent in mydict[spectralSpec]:
                spectralSpec += 1
                if spectralSpec not in mydict:
                    mydict[spectralSpec] = {}
            if intent not in list(mydict[spectralSpec].keys()):
                mydict[spectralSpec][intent] = {'sourcename': sourcename, 'originalIntegration': float(originalIntegration), 
                                                'optimumIntegration': float(optimumIntegration), 'adoptedIntegration': float(adoptedIntegration)}
    return mydict

def readSoftwareVersionFromASDM_minidom(asdm):
    """
    Reads the software version from the ASDM's Annotation.xml table.
    - Todd Hunter
    """
    if (os.path.exists(asdm) == False):
        print("readSoftwareVersionFromASDM_minidom(): Could not find ASDM = ", asdm)
        return(None)
    if (os.path.exists(asdm+'/Annotation.xml') == False):
        print("readSoftwareVersionFromASDM_minidom(): Could not find Annotation.xml. This dataset was probably taken prior to R10.6.")
        return(None)

    xmlscans = minidom.parse(asdm+'/Annotation.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    print('\n### Software version for ASDM: %s ###' % asdm)
    for i,rownode in enumerate(rowlist):
        row = rownode.getElementsByTagName("issue")
        issue = str(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("details")
        details = str(row[0].childNodes[0].nodeValue)
        print("%s: %s" % (issue,details))
    return

def getWVREfficienciesFromASDM(asdm):
    mydict = {}
    xml = asdm + '/CalReduction.xml'
    xmlscans = minidom.parse(xml)
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    mydict = {}
    for rownode in rowlist:
        row = rownode.getElementsByTagName("paramSet")
        parameterSet = str(row[0].childNodes[0].nodeValue)
        if (parameterSet.find('serialNumber') >= 0):
            phrases = parameterSet.split()[2:]
            for phrase in phrases:
                serialNumber = phrase.find('serialNumber')
                skyCoupling = phrase.find('skyCoupling')
                if (serialNumber >= 0):
                    antenna = phrase[serialNumber+13:serialNumber+17]
                    if antenna not in mydict.keys():
                        mydict[antenna] = {}
                    mydict[antenna]['serialNumber'] = int(phrase[serialNumber+19:].rstrip('"'))
                elif (skyCoupling >= 0):
                    antenna = phrase[skyCoupling+12:skyCoupling+16]
                    if antenna not in mydict.keys():
                        mydict[antenna] = {}
                    mydict[antenna]['skyCoupling'] = float(phrase[skyCoupling+19:].rstrip('"'))
                    
    return mydict

def getAntennaPadsFromASDM_minidom(asdm):
    mydict = readStationFromASDM_minidom(asdm)
    pads = []
    for key in mydict.keys():
        pad = mydict[key]['name'].strip() 
        padtype = mydict[key]['type']
        if (padtype == ANTENNA_PAD):
            pads.append(pad)
#        if (pad.find('WSTB') < 0):
#            pads.append(pad)
    return pads

def readAntennaPositionFromASDM_minidom(sdmfile, antennaType=''):
    """
    Reads the Antenna.xml file and returns a dictionary of all antennas
    of the following format:
    mydict = {'DV04': {'id': 0, 'position': [x,y,z]}}
    -Todd Hunter
    """
    if (os.path.exists(sdmfile) == False):
        print("readAntennaPositionFromASDM_minidom(): Could not find file = ", sdmfile)
        return(None)
    xmlscans = minidom.parse(sdmfile+'/Antenna.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    stationName = 'unknown'
    mydict = {}
    positions = []
    for rownode in rowlist:
        stationPosition = []
        scandict[fid] = {}
        row = rownode.getElementsByTagName("antennaId")
        stationId = int(str(row[0].childNodes[0].nodeValue).split('_')[-1])
        row = rownode.getElementsByTagName("name")
        stationName = str(row[0].childNodes[0].nodeValue).strip()
        row = rownode.getElementsByTagName("position")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        for i in range(2,len(r)):
            stationPosition.append(float(r[i]))
        if antennaType == '' or stationName.find(antennaType)==0:
            mydict[stationName] = {'id': fid, 'position': stationPosition}
            fid +=1
            positions.append(stationPosition)
    if antennaType != '':
        positions = np.array(positions)
        medianVector = np.median(positions, axis=0)
        positions = np.transpose(positions-medianVector)
        print("median position: X=%+f Y=%+f Z=%+f" % (medianVector[0],medianVector[1],medianVector[2]))
        print("rms variation:   X=%+f Y=%+f Z=%+f" % (np.std(positions[0]),np.std(positions[1]),np.std(positions[2])))
    return(mydict)
    
def readStationFromASDM_minidom(sdmfile):
    """
    Reads the Station.xml file and returns a dictionary of all stations
    of the following format:
    mydict[0] = {'name': 'A085', 'position': [x,y,z]}
    -Todd Hunter
    """
    if (os.path.exists(sdmfile) == False):
        print("readStationFromASDM_minidom(): Could not find file = ", sdmfile)
        return(None)
    xmlscans = minidom.parse(sdmfile+'/Station.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    stationName = 'unknown'
    mydict = {}
    for rownode in rowlist:
        stationPosition = []
        scandict[fid] = {}
        row = rownode.getElementsByTagName("stationId")
        stationId = int(str(row[0].childNodes[0].nodeValue).split('_')[-1])
        row = rownode.getElementsByTagName("name")
        stationName = str(row[0].childNodes[0].nodeValue).strip()
        row = rownode.getElementsByTagName("type")
        stationType = str(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("position")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        for i in range(2,len(r)):
            stationPosition.append(float(r[i]))
        mydict[stationId] = {'name': stationName, 'position': stationPosition, 'type': stationType}
        fid +=1
    return(mydict)
    
def readStationsFromASDM_minidom(sdmfile, station=None):
    """
    Translates a station number (which start from 0) into the station name and
    position from the Station.xml file.  Useful for finding this information
    for weather stations.
    If station==None, then it builds and returns a dictionary where the key is
    the station name and the value is the geocentric [X,Y,Z] position.
    e.g. {'A001': [x,y,z]}
    - Todd Hunter
    """
    if (os.path.exists(sdmfile) == False):
        print("readStationFromASDM()_minidom: Could not find file = ", sdmfile)
        return(None)
    xmlscans = minidom.parse(sdmfile+'/Station.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    fid = 0
    stationName = 'unknown'
    if (station == None):
        mydict = {}
    for rownode in rowlist:
        stationPosition = []
        scandict[fid] = {}
        row = rownode.getElementsByTagName("stationId")
        stationId = int(str(row[0].childNodes[0].nodeValue).split('_')[-1])
        row = rownode.getElementsByTagName("name")
        stationName = str(row[0].childNodes[0].nodeValue).strip() # remove spaces added in 2017 by Oracle
        row = rownode.getElementsByTagName("position")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        for i in range(2,len(r)):
            stationPosition.append(float(r[i]))
        if (stationId == station):
            break
        elif (station == None):
            mydict[stationName] = stationPosition
        fid +=1
    if (station == None):
        return(mydict)
    else:
        return(stationName,stationPosition)

def getSubscanTimesFromASDM_minidom(asdm, field=''):
    """
    Reads the subscan information from the ASDM's Subscan.xml file and
    returns a dictionary of form:
    {scan: {subscan: {'field': '3c273, 'integrationTime': 2.016,
                      'numIntegration': 5, 'subscanLength': 10.08}}}
    where the scan numbers are the top-level keys.  The subscanLength is
    computed by the difference between endTime and startTime.  The integration
    time is computed by dividing the subscanLength by numIntegration.
    If the field name is specified, then limit the output to scans on this
    field.
    -- Todd Hunter
    """
    subscanxml = asdm + '/Subscan.xml'
    if (os.path.exists(subscanxml) == False):
        print("Could not open %s" % (subscanxml))
        return
    xmlscans = minidom.parse(subscanxml)
    rowlist = xmlscans.getElementsByTagName("row")
    scandict = {}
    scanNumbers = 0
    subscanTotalLength = 0
    for rownode in rowlist:
        row = rownode.getElementsByTagName("scanNumber")
        scanNumber = int(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("subscanNumber")
        subscanNumber = int(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("startTime")
        startTime = int(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("endTime")
        endTime = int(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("numIntegration")
        numIntegration = int(row[0].childNodes[0].nodeValue)
        row = rownode.getElementsByTagName("fieldName")
        fieldName = str(row[0].childNodes[0].nodeValue)
        if (field=='' or fieldName==field):
            subscanLength = (endTime-startTime)*1e-9
            subscanTotalLength += subscanLength
            integrationTime = subscanLength / (1.0*numIntegration)
            if (scanNumber not in scandict):
                if (scanNumber == 1):
                    scan1startTime = startTime
                scandict[scanNumber] = {}
                scanNumbers += 1
            scandict[scanNumber][subscanNumber] = {'subscanLength': subscanLength, 'numIntegration': numIntegration, 'integrationTime': integrationTime, 'field': fieldName, 'startTime':startTime*1e-9, 'endTime':endTime*1e-9}
    print("Found %d scans" % (scanNumbers))
    totalTime = (endTime-scan1startTime)*1e-9
    latency = totalTime - subscanTotalLength
    print("Total latency = %g/%g seconds = %g percent" % (latency, totalTime, latency*100/totalTime))
    return(scandict)

def readDecorrelationFromASDM_minidom(asdm):
    """
    -Todd Hunter
    """
    mydict = {}
    seeingxml = asdm + '/CalPhase.xml'
    if (os.path.exists(seeingxml) == False):
        print("Could not open %s" % (seeingxml))
        return
    xml = minidom.parse(seeingxml)
    rowlist = xml.getElementsByTagName("row")
    mydict['basebandName'] = []
    mydict['receiverBand'] = []
    mydict['numReceptor'] = []
    mydict['baselineLengths'] = []
    mydict['decorrelationFactor'] = []
    mydict['startValidTime'] = []
    mydict['endValidTime'] = []
    mydict['atmPhaseCorrection'] = []
    mydict['integrationTime'] = []
    mydict['azimuth'] = []
    mydict['elevation'] = []
    mydict['calDataId'] = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("startValidTime")
        mydict['startValidTime'].append(int(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("endValidTime")
        mydict['endValidTime'].append(int(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("atmPhaseCorrection")
        mydict['atmPhaseCorrection'].append(str(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("receiverBand")
        mydict['receiverBand'].append(str(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("basebandName")
        mydict['basebandName'].append(str(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("numReceptor")
        numReceptor = int(row[0].childNodes[0].nodeValue)
        mydict['numReceptor'].append(numReceptor)
        row = rownode.getElementsByTagName("calDataId")
        mydict['calDataId'].append(int(str(row[0].childNodes[0].nodeValue).split('_')[1]))
        row = rownode.getElementsByTagName("integrationTime")
        mydict['integrationTime'].append(float(row[0].childNodes[0].nodeValue)*1e-9)
        row = rownode.getElementsByTagName("baselineLengths")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        baselineLengths = []
        for i in range(2,len(r)):
            baselineLengths.append(float(r[i]))
        mydict['baselineLengths'].append(baselineLengths)
        row = rownode.getElementsByTagName("decorrelationFactor")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        decorrelationFactor = []
        for i in range(3,len(r)):
            decorrelationFactor.append(float(r[i]))
        mydict['decorrelationFactor'].append(decorrelationFactor)
        row = rownode.getElementsByTagName("direction")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        direction = []
        for i in range(2,len(r)):
            direction.append(float(r[i]))
        mydict['azimuth'].append(math.degrees(direction[0]))
        mydict['elevation'].append(math.degrees(direction[1]))
    print("Found %d measurements on %d baselines" % (len(mydict['atmPhaseCorrection']), len(mydict['baselineLengths'][0])))
    return mydict
    
def readSeeingFromASDM_minidom(asdm):
    """
    Reads information from CalSeeing.xml into a dictionary
    Returns a dictionary with the following keys:
    atmPhaseCorrection: AP_UNCORRECTED or AP_CORRECTED
    baselineLengths: typically 3 values (in meters)
    startValidTime: MJD nano seconds
    endValidTime: MJD nano seconds
    phaseRMS:  a value for each baselineLength (radians?) for each timestamp
    seeing: one value per timestamp (arcseconds)
    -Todd Hunter
    """
    mydict = {}
    seeingxml = asdm + '/CalSeeing.xml'
    if (os.path.exists(seeingxml) == False):
        print("Could not open %s" % (seeingxml))
        return
    xml = minidom.parse(seeingxml)
    rowlist = xml.getElementsByTagName("row")
    mydict['seeing'] = []
    mydict['phaseRMS'] = []
    mydict['startValidTime'] = []
    mydict['endValidTime'] = []
    mydict['atmPhaseCorrection'] = []
    mydict['baselineLengths'] = []
    mydict['phaseRMS'] = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("seeing")
        mydict['seeing'].append(float(row[0].childNodes[0].nodeValue)*206264.8)
        row = rownode.getElementsByTagName("startValidTime")
        mydict['startValidTime'].append(int(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("endValidTime")
        mydict['endValidTime'].append(int(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("atmPhaseCorrection")
        mydict['atmPhaseCorrection'].append(str(row[0].childNodes[0].nodeValue))
        row = rownode.getElementsByTagName("baselineLengths")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        baselineLengths = []
        for i in range(2,len(r)):
            baselineLengths.append(float(r[i]))
        mydict['baselineLengths'].append(baselineLengths)
        row = rownode.getElementsByTagName("phaseRMS")
        r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
        phaseRMS = []
        for i in range(2,len(r)):
            phaseRMS.append(float(r[i]))
        mydict['phaseRMS'].append(phaseRMS)
    print("Found %d measurements" % (len(mydict['atmPhaseCorrection'])))
    return mydict

def asdmspwmap_minidom(asdm):
    """
    Generate a list that maps the spw number that will be found in the
    measurement set to the corresponding value in the ASDM xml files.
    In general, the order will be [0,n+1,n+2,....] where n=number of antennas
    with WVR data.  For example: [0,5,6,7...] if n=4 antennas, meaning
    that spw 1 in the ms = spw 5 in the ASDM xml files.
    -Todd Hunter
    """
    mydict = readSpwsFromASDM_minidom(asdm)
    spws = []
    for i,spw in enumerate(mydict['spw']):
        if (mydict['name'][i].find('WVR#Antenna') < 0):
            spws.append(int(i))
    return(spws)

def readSpwsFromASDM_minidom(asdm, verbose=False):
    """
    Reads spw information from SpectralWindow.xml into a dictionary
    Returns a dictionary with the following keys:
    'spw': string number
    'name': string e.g. 'WVR#NOMINAL'
    -Todd Hunter
    """
    mydict = {}
    wvrAntennas = 0
    spwxml = asdm + '/SpectralWindow.xml'
    if (os.path.exists(spwxml) == False):
        print("Could not open %s" % (spwxml))
        return
    xml = minidom.parse(spwxml)
    rowlist = xml.getElementsByTagName("row")
    mydict['spw'] = []
    mydict['name'] = []
    for rownode in rowlist:
        row = rownode.getElementsByTagName("name")
        name = str(row[0].childNodes[0].nodeValue)
        mydict['name'].append(name)
        row = rownode.getElementsByTagName("spectralWindowId")
        mydict['spw'].append(str(row[0].childNodes[0].nodeValue).split('_')[1])
        if (name.find('#Antenna') > 0):
            wvrAntennas += 1
    if verbose:
        print("Found %d spws" % (len(mydict['spw'])))
        if (wvrAntennas > 0):
            print("but %d are only for the WVR filter frequencies." % (wvrAntennas))
    return mydict

def readFluxesFromASDM_minidom(sdmfile, useCalFlux=False, sourcename='', field=-1, spw=-1):
    if (useCalFlux):
        calflux = sdmfile + '/CalFlux.xml'
        if not os.path.exists(calflux):
            print("readFluxesFromASDM_minidom(): Could not find file = ", calflux)
            print("Looking for Source.xml instead")
            sourcefile = sdmfile + '/Source.xml'
            if (os.path.exists(sourcefile) == False):
                print("readFluxesFromASDM_minidom(): Could not find file = ", sourcefile)
                return(None)
            useCalFlux = False
    else:
        sourcefile = sdmfile + '/Source.xml'
        if not os.path.exists(sourcefile):
            print("readFluxesFromASDM(): Could not find file = ", sourcefile)
            return(None)
    scandict = {}
    fid = 0
    sources = []
    myspw = ''
    sourceid=''
    noFrequencyKeywords = True
    if not useCalFlux:   # Source.xml
        sourcesIdentified = {}
        xmlscans = minidom.parse(sourcefile)
        rowlist = xmlscans.getElementsByTagName("row")
        for rownode in rowlist:
            scandict[fid] = {}
            row = rownode.getElementsByTagName("sourceId")
            sourceId = int(str(row[0].childNodes[0].nodeValue))
            if (sourceid != '' and sourceid != sourceId): continue
            row = rownode.getElementsByTagName("sourceName")
            sourceName = str(row[0].childNodes[0].nodeValue).strip(' ')
            if sourceId not in sourcesIdentified:
                sourcesIdentified[sourceId] = [sourceName]
            elif sourceName not in sourcesIdentified[sourceId]:
                print("WARNING: There is a mismatch in the sourceID/sourceName. ID%d has multiple names: " % (sourceId), sourceName, sourcesIdentified[sourceId])
                sourcesIdentified[sourceId].append(sourceName)
            if (sourcename != ''):
                if (sourcename.find('*') >= 0):
                    if (sourceName.find(sourcename.replace('*','')) < 0): continue
                else:
                    if (sourceName != sourcename):
                        continue
            elif (int(field) != -1):
                if (sourceId != int(field)):
                    continue
                
            row = rownode.getElementsByTagName("spectralWindowId")
            spectralWindowId = int(str(row[0].childNodes[0].nodeValue).split('_')[-1])
            if (myspw != '' and myspw != spectralWindowId): continue
            row = rownode.getElementsByTagName("frequency")
            if (row == []):
                continue
            noFrequencyKeywords = False
            r = list(filter(None, (row[0].childNodes[0].nodeValue).split(' ')))
            frequency = []
            for freq in r[2:]:
                frequency.append(float(freq))
            row = rownode.getElementsByTagName("stokesParameter")
            r = list(filter(None,str(row[0].childNodes[0].nodeValue).split(' ')))
            stokesParameters = []
            for stokes in r[2:]:
                stokesParameters.append(stokes)
            row = rownode.getElementsByTagName("flux")
            r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
            nflux = int(r[1])
            fluxDensity = []
            for flux in range(nflux):
                stokesFlux = []
                for stokes in range(len(stokesParameters)):
                    stokesFlux.append(float(r[3+stokes+flux*len(stokesParameters)]))
                fluxDensity.append(stokesFlux)
            if (spectralWindowId in spw or spw==[-1]):
                sources.append({'sourceId': sourceId, 'sourceName': sourceName,
                                'spw': spectralWindowId,
                                'stokesParameters': stokesParameters,
                                'frequency': frequency, 'fluxDensity': fluxDensity})
                fid +=1
    else:  # CalFlux.xml
        xmlscans = minidom.parse(calflux)
        rowlist = xmlscans.getElementsByTagName("row")
        for rownode in rowlist:
            scandict[fid] = {}
            row = rownode.getElementsByTagName("sourceName")
            sourceName = str(row[0].childNodes[0].nodeValue)
            if (sourcename != ''):
                if (sourcename.find('*') >= 0):
                    if (sourceName.find(sourcename.replace('*','')) < 0): continue
                else:
                    if (sourceName != sourcename):
                        continue
            row = rownode.getElementsByTagName("frequencyRanges")
            if (row == []):
                continue
            r = list(filter(None, (row[0].childNodes[0].nodeValue).split(' ')))
            frequency = []
            for freq in range(3,len(r)-1,2):
                frequency.append([float(r[freq])*1.0e-9, float(r[freq+1])*1.0e-9])
            row = rownode.getElementsByTagName("stokes")
            r = list(filter(None,str(row[0].childNodes[0].nodeValue).split(' ')))
            stokesParameters = []
            for stokes in r[2:]:
                stokesParameters.append(stokes)
            row = rownode.getElementsByTagName("flux")
            r = list(filter(None,(row[0].childNodes[0].nodeValue).split(' ')))
            rowerror = rownode.getElementsByTagName("fluxError")
            rerror = list(filter(None,(rowerror[0].childNodes[0].nodeValue).split(' ')))
            nflux = int(r[1])
            print("nflux = ", nflux)
            fluxDensity = []
            fluxDensityError = []
            for flux in range(nflux):
                stokesFlux = []
                errorFlux = []
                for stokes in range(len(stokesParameters)):
                    stokesFlux.append(float(r[3+stokes+flux*len(stokesParameters)]))
                    errorFlux.append(float(rerror[3+stokes+flux*len(stokesParameters)]))
                fluxDensity.append(stokesFlux)
                fluxDensityError.append(errorFlux)
            sources.append({'sourceName': sourceName, 'stokes': stokesParameters,
                            'frequencyRange': frequency, 'fluxDensity': fluxDensity,
                            'fluxDensityError':fluxDensityError})
            fid +=1
    return sources, noFrequencyKeywords, useCalFlux

def createCasaTool(mytool):
    """
    A wrapper to handle the changing ways in which casa tools are invoked.
    For CASA < 6, it relies on "from taskinit import *" in the preamble above.
    mytool: a tool name, like tbtool
    Todd Hunter
    """
    if 'casac' in mylocals:
        if (type(casac.Quantity) != type):  # casa 4.x and 5.x
            myt = mytool()
        else:  # casa 3.x
            myt = mytool.create()
    else:
        # this is CASA 6
        myt = mytool()
    return(myt)

