#!/usr/bin/python
"""
Script to run the obspy picking routines, takes a few arguments of the form:
Starttime, endtime, dataloc, routine type (classic, carl)
Data must be formatted in hour long multiplexed miniseed files and named in
seisan yyyy-mm-dd-hhmm-ss format

v 1.0 - Working, but parameters set to defaults
v 1.1 - Outputs to seisan directory structure
Calum Chamberlain, Victoria University of Wellington, 09/11/2014
"""
def cjc_trigger_routine(startdate,enddate,dataloc,trigloc,routype):
    """
    Module to run the obspy sta-lta energy based filter routine

    Must be parsed start date & end date in obspy UTCDateTime type,
    dataloc should be a string of the path for the input data
    trigloc should be a string of the ouput path
    routype should be a string denpoting the type of detection routine to use
        either classic or carl
    defaults have been set in the module for trigger parameters
    """

###############################################################################
    # Import parameter settings
    import sys
    sys.path.insert(0,"/home/calumch/my_programs/Building/rt2detection")
    from par import trigger_par as defaults
    print defaults.stalen
###############################################################################

# Format dates
    startyear=startdate.split('/')[0]
    startmonth=startdate.split('/')[1]
    startday=startdate.split('/')[2]
    endyear=enddate.split('/')[0]
    endmonth=enddate.split('/')[1]
    endday=enddate.split('/')[2]

# Import modules
    from obspy import read as obsread
    from obspy import UTCDateTime
    import glob, os
    import numpy as np
    from obspy.signal import coincidenceTrigger


# Generate list of days to check through
    lengthinseconds=UTCDateTime(endyear+' '+endmonth+' '+endday)-\
            UTCDateTime(startyear+' '+startmonth+' '+startday)
    lendays=lengthinseconds/86400
    lengthinseconds=[]
    dfiles=[]
    dates=[]
    for i in range(0,int(lendays)+1):
        dates.append(UTCDateTime(startyear+' '+startmonth+' '+startday)+(i*86400))
        dfiles.extend(glob.glob(dataloc+'/'+str(dates[i].year)+'/'+\
                str(dates[i].month).zfill(2)+'/'+str(dates[i].year)+'-'+\
                str(dates[i].month).zfill(2)+'-'+str(dates[i].day).zfill(2)+'*'))

    print len(dfiles)
    wavelist=[] # Initialize list variable
    # Read in data
    for hfile in dfiles:
        print 'Working on file: '+hfile
        st=obsread(hfile)
        st1=st.copy()
        if not defaults.comp=='all':
            st1=st1.select(channel='*'+defaults.comp)
        # De-mean data
        for tr in st:
            tr.data=tr.data-np.mean(tr.data)
        # Filter data
        st1.filter('bandpass',freqmin=defaults.lowcut,freqmax=defaults.highcut)
        # Use the obspy triggering routine
        trig=[]
        if routype=='classic':
            trig = coincidenceTrigger("recstalta",defaults.trigon,\
                                      defaults.trigoff,st1,defaults.netsum,\
                                      sta=defaults.stalen,lta=defaults.ltalen,\
                                      delete_long_trigger='True',\
                                      trigger_off_extension=\
                                      defaults.netwin)
        else:
            try:
                trig = coincidenceTrigger("carlstatrig",defaults.trigon,\
                                          defaults.trigoff,st1,\
                                          defaults.netsum,sta=defaults.stalen,\
                                          lta=defaults.ltalen,ratio=defaults.crat,\
                                          quiet=defaults.cquite,delete_long_trigger='True')
            except:
                print 'Triggering routine failed, suggest altering parameters'
        # Cut data and write out in multiplexed miniseed files
        if trig and defaults.trigout=='Y':
            for event in trig:
                stout=st.slice(event['time']-defaults.precut,event['time']+defaults.postcut)
                filename=str(stout[0].stats.starttime.year)+'-'+\
                        str(stout[0].stats.starttime.month).zfill(2)+'-'+\
                        str(stout[0].stats.starttime.day).zfill(2)+'-'+\
                        str(stout[0].stats.starttime.hour).zfill(2)+\
                        str(stout[0].stats.starttime.minute).zfill(2)+'-'+\
                        str(stout[0].stats.starttime.second).zfill(2)+'.'+\
                        defaults.net+'_'+str(len(stout)).zfill(3)+'_00'
                if not os.path.isdir(trigloc+'/'+\
                        str(stout[0].stats.starttime.year)):
                    os.makedirs(trigloc+'/'+str(stout[0].stats.starttime.year))
                if not os.path.isdir(trigloc+'/'+str(stout[0].stats.starttime.year)\
                        +'/'+str(stout[0].stats.starttime.month).zfill(2)):
                    os.makedirs(trigloc+'/'+str(stout[0].stats.starttime.year)\
                            +'/'+str(stout[0].stats.starttime.month).zfill(2))
                filename=trigloc+'/'+str(stout[0].stats.starttime.year)+'/'+\
                        str(stout[0].stats.starttime.month).zfill(2)+'/'+\
                        filename
                wavelist.append(filename)
                try:
                    stout.write(filename,format="MSEED",encoding="STEIM2")
                except:
                    # Cope with dtype issues
                    for tr in stout:
                        tr.data = np.array(tr.data, dtype=np.int32)
                    stout.write(filename,format='MSEED',encoding='STEIM2')
                print 'Written triggered file as: '+filename
        elif defaults.trigout=='N':
            print 'Triggers will not be written out but I made '+len(trig)+' detections'
        elif not trig:
            print 'No triggers were detected'
    return wavelist

if __name__=='__main__':
    # Read arguments and check
    import sys, os
    if len(sys.argv) != 6:
        print 'Requires 5 arguments of: start date (yyyy/mm/dd), end date (yyyy/mm/dd)'
        print 'data location in seisan style folders, output location in seisan'
        print ' style folders, and routine type (classic or carl)'
        sys.exit()
    else:
        startdate=str(sys.argv[1])
        if len(startdate.split('/')) != 3:
            print 'Startdate must be formatted yyyy/mm/dd'
            sys.exit()
        enddate=str(sys.argv[2])
        if len(enddate.split('/')) != 3:
            print 'End date must be formatted yyyy/mm/dd'
            sys.exit()
        dataloc=str(sys.argv[3])
        if not os.path.isdir(dataloc):
            print 'Path given for seismic data does not exist'
            print dataloc
            sys.exit()
        trigloc=sys.argv[4]
        if not os.path.isdir(trigloc):
            print 'Path given for output does not exist'
            print trigloc
            sys.exit()
        routype=str(sys.argv[5])
        if routype not in ['classic','carl']:
            print 'Must give either classic or carl as final argument'
            sys.exit()

    wavelist=cjc_trigger_routine(startdate,enddate,dataloc,trigloc,routype)
    if wavelist==[]:
        print 'No detections made'
    else:
        print 'Triggered file list:'
    for wavfile in wavelist:
        print wavfile
