#!/usr/bin/python
"""
Code designed to take data downloded using defaults.neo and located in the RT_data
folder of this ditribution and convert this to miniseed, multiplex it and run
seisan's continuous detection, sta/lta routine over this data.

ver 0.1 - Converts and defaults.merges data, but only for one day

ver 0.2 - Converts and defaults.merges data for multiple days

ver 1.1 - Converts and defaults.merges data for multiple days and can run a python
         detection routine (STA/LTA) over the data

ver 2.1 - Uses the filterdefaults.picker routine to automatically pick data, parameters
         for this are not yet optimized

07/11/14 - Calum Chamberlain - Victoria University of Wellington, NZ
"""
############################ FUNCTIONS ########################################

def plot_data_cont(stations, outfile, form):
    import datetime as dt
    import matplotlib.dates as mpdates
    import matplotlib.pyplot as plt
    plotno=1 # counter to move plot line up
    fig=plt.figure(num=1, dpi=100, facecolor='w', edgecolor='k')
    fig.suptitle('Data Continuity')
    # First work out how long the network has been active for
    mindate=mpdates.date2num(dt.datetime.strptime('2100365','%Y%j'))
    maxdate=mpdates.date2num(dt.datetime.strptime('1900365','%Y%j'))
    for sta in stations:
        f = open(sta.name+'_daslist.csv','r')
        date=[]
        for line in f:
            date+=[line.split(',')[0]]
        f.close()
        dates=[]
        for d in date:
            d=d[0:7]
            dates+=[dt.datetime.strptime(d,'%Y%j')]
        dates=sorted(list(set(dates)))
        dates=mpdates.date2num(dates)
        if dates[0] < mindate:
            mindate=dates[0]
        if dates[len(dates)-1] > maxdate:
            maxdate=dates[len(dates)-1]
    # Generate X axis values - consistent throughout
    x=range(int(mindate),int(maxdate))
    xdates=[mpdates.num2date(xd) for xd in x]
    # Now re-read in data for useful-ness: almost certainly not the most efficient
    # way to do this
    for sta in stations:
        f = open(sta.name+'_daslist.csv','r')
        date=[]
        for line in f:
            date+=[line.split(',')[0]]
        f.close()
        dates=[]
        for d in date:
            d=d[0:7]
            dates+=[dt.datetime.strptime(d,'%Y%j')]
        dates=sorted(list(set(dates)))
        dates=mpdates.date2num(dates)
        # convert dates to int
        dates=[int(d) for d in dates]
        y=[]
        for d in x:
            if d in dates:
                y+=[1]
            else:
                y+=[0]
        # Here is where the plotting happens
        ax = fig.add_subplot(len(stations), 1, plotno)
        ax.scatter(xdates,y)
        ax.set_ylim([0.5,1.5])
        ax.set_xlim([min(x),max(x)])
        ax.set_ylabel(sta.name, rotation=0)
        ax.xaxis.set_major_locator(mpdates.YearLocator())
        ax.xaxis.set_major_formatter(mpdates.DateFormatter('%Y-%m-%d'))
        ax.xaxis.set_minor_locator(mpdates.MonthLocator())
        plotno+=1
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
    ax.fmt_xdata = mpdates.DateFormatter('%Y-%m-%d')
    fig.autofmt_xdate()
    fig.set_size_inches(12,10)
    plt.savefig(outfile, format=form)

if __name__ == '__main__':
    ## GET DEFAULT PARAMETERS
    from par import overallpars as defaults
    if defaults.debug==1:
        print "Looking for data in: "+defaults.indir
        print "Will put miniseed data in: "+defaults.outdir
        print "Planning to put merged data here: "+defaults.contbase
        print "Planning to put s-files here: "+defaults.sfilebase
        print "Planning to put triggered waveforms here: "+defaults.trigbase
        print "The network code you have selected is: "+defaults.network
        print "I will use the following detection routine: "+defaults.detrout
        print "Which will run the following style of detection: "+defaults.routype
        print "Continuous waveforms will be of this length (s): "+str(defaults.filelenwant)
        print "Continuous waveforms will be of this mseed format: "+defaults.dformat
        if not defaults.rmold:
                print "I will not remove old data, may result in me doing more than I need to"
        else:
            print " I will be removing old data, watch yourself!"
        print "I will pick the data using the following routine: "+defaults.picker
        if defaults.rawconv:
            print "I will be converting data from raw today"
        else:
            print "No raw conversion selected, data must be in ms format in the correct folder"
        if defaults.merge:
            print "Data will be merged into multiplexed miniseed files"
            print "Maximum residual allowed for picker: "+str(defaults.maxres)
            print "User ID: "+defaults.userID
        if not defaults.overwrite:
            print "Will conserve old s-files"
        else:
            print "Will overwrite old s-files"
        if defaults.neo=='True':
            print "Expecting neo .ZIP files"
        if not defaults.alltime:
            print 'Will not run for all the data'
            print 'Startime is: '+defaults.startD+', end time is: '+defaults.endD
        print "\n--------------------------------------------------------------\n\n"
        if not defaults.alltime:
            from obspy import UTCDateTime
            startD=UTCDateTime(int(defaults.startD[0:4]),int(defaults.startD[5:7]),\
                               int(defaults.startD[8:10]))
            endD=UTCDateTime(int(defaults.endD[0:4]),int(defaults.endD[5:7]),\
                             int(defaults.endD[8:10]))
            rundates=[]
            for i in range(0,int((endD-startD)/86400)):
                rundates+=[startD+(i*86400)]
        from par.overallpars import STATION

    # Import necessary libraries
    from subprocess import call
    import sys, datetime
    from copy import deepcopy

    # Open log-file
    f = open('rt2detection_log.txt','w')

        ############################CODE START#########################################

    # Clear old miniseed data so that we don't read in a ridiculous amount of data
    import shutil, glob, os
    if defaults.rawconv and defaults.rmold:
        rmlist=glob.glob(defaults.outdir+'/*')
        for tree in rmlist:
            print 'Removing old data'
            shutil.rmtree(tree)

    ################### Convert Reftek to miniseed#################################
    # This section of the code requires you to have archived the raw data using
    # pascal tools defaults.neo tool in station directories - it then uses the pascal tool
    # rt2ms to convert from these archived data to miniseed files in the length
    # that they were recorded in.
    if defaults.rawconv and defaults.neo=='True':
        for sta in defaults.stalist:
            print('Running rt2ms for station: '+sta.name)
            call(["rt2ms","-n",sta.netcode,"-s",sta.name,"-l",sta.loccode,\
                    "-d",defaults.indir+"/"+sta.name,"--daskeep",sta.das,"-o",\
                    defaults.outdir+'/'+sta.name,"-Y"])
    if defaults.rawconv and defaults.neo=='False':
        for sta in defaults.stalist:
            print('Generating lists of the files to be converted for station '+\
                  sta.name)
            if not defaults.alltime:
                filelist=[]
                daslist=[]
                for date in rundates:
                    datestr=str(date.year)+str(date.julday).zfill(3)
                    filelist+=(glob.glob(defaults.indir+'/'+sta.name+'/'+datestr))
                    daslist+=(glob.glob(defaults.indir+'/'+sta.name+'/'+datestr+'/*'))
                filelist.sort()
                daslist.sort()
            else:
                filelist=sorted(glob.glob(defaults.indir+'/'+sta.name+'/*'))
                daslist=sorted(glob.glob(defaults.indir+'/'+sta.name+'/*/*'))
            if len(filelist) != 0:
                rtflist = open('rtflist.lis','w')
                for rtfile in filelist:
                    rtflist.write(rtfile+'\n')
                rtflist.close()
                dasfilelist = open(sta.name+'_daslist.csv','w')
                for das in daslist:
                    dasfilelist.write(das.split('/')[len(das.split('/'))-2]+','+\
                                        das.split('/')[len(das.split('/'))-1]+'\n')
                dasfilelist.close()
                call(["rt2ms","-n",sta.netcode,"-s",sta.name,"-l",sta.loccode,\
                        "-F","rtflist.lis","--daskeep",sta.das,"-o",\
                        defaults.outdir+'/'+sta.name,"-Y","-L"])
            else:
                print('Found no data for station: '+sta.name)

    if defaults.rawconv and defaults.neo=='Sometimes':
        # Special case to look for both files
        for sta in defaults.stalist:
            print('Running rt2ms on neo input files first for station: '+sta.name)
            # call(["rt2ms","-n",sta.netcode,"-s",sta.name,"-l",sta.loccode,\
                    # "-d",defaults.indir+"/"+sta.name,"--daskeep",sta.das,"-o",\
                    # defaults.outdir+'/'+sta.name,"-Y"])
            print('Generating list of files to be run through rt2ms for station: '+\
                  sta.name)
            if not defaults.alltime:
                filelist=[]
                daslist=[]
                for date in rundates:
                    filelist+=(glob.glob(defaults.indir+'/'+sta.name+'/'+date))
                    daslist+=(glob.glob(defaults.indir+'/'+sta.name+'/'+date+'/*'))
                filelist.sort()
                daslist.sort()
            else:
                filelist=sorted(glob.glob(defaults.indir+'/'+sta.name+'/*'))
                daslist=sorted(glob.glob(defaults.indir+'/'+sta.name+'/*/*'))
            if len(filelist) != 0:
                rtflist = open('rtflist.lis','w')
                for rtfile in filelist:
                    if not rtfile[len(rtfile)-4:len(rtfile)] == '.ZIP':
                        rtflist.write(rtfile+'\n')
                    else:
                        daslist+=rtfile
                rtflist.close()
                dasfilelist = open(sta.name+'_daslist.csv','w')
                for das in daslist:
                    dasfilelist.write(das.split('/')[len(das.split('/'))-2]+','+\
                                        das.split('/')[len(das.split('/'))-1]+'\n')
                dasfilelist.close()
                # call(["rt2ms","-n",sta.netcode,"-s",sta.name,"-l",sta.loccode,\
                        # "-F","rtflist.lis","--daskeep",sta.das,"-o",\
                        # defaults.outdir+'/'+sta.name,"-Y","-L"])
            else:
                print('Found no data for station: '+sta.name)

    ############################ Plot data continuity #############################
    if defaults.rawconv or defaults.archive:
        try:
            plot_data_cont(defaults.stalist,'Data_continuity.eps','eps')
        except:
            print 'Did not find station continuity files'

    ######################## Archive data #########################################
    from obspy import read as obsread
    from obspy import UTCDateTime, Stream, Trace
    ######################## merge data ###########################################
    # merge miniseed data to hour long files, this requires obspy to easily read
    # and manipulate the data.

    # defaults.merge data to a multiplexed files

    if defaults.merge or defaults.archive:
        # Check that output does not span multiple days
        if not defaults.alltime:
            if defaults.rawconv:
                daylist=[]
                for date in rundates:
                    daylist+=glob.glob(defaults.outdir+'/*/Y'+str(date.year)+\
                                       '/R'+str(date.julday).zfill(3))+'.01'
            else:
                daylist=[]
                for date in rundates:
                    daylist+=glob.glob(defaults.outdir+'/Y'+str(date.year)\
                                      +'/R'+str(date.julday).zfill(3)+'.01')
        else:
            if defaults.rawconv:
                daylist=glob.glob(defaults.outdir+'/*/Y*/R*')
            else:
                daylist=glob.glob(defaults.outdir+'/Y*/R*')
        if defaults.debug==1:
            print 'Found '+str(len(daylist))+' days to run through'
        datelist=[daylist[0].split('/')[len(daylist[0].split('/'))-2]+'/'+\
                    daylist[0].split('/')[len(daylist[0].split('/'))-1]]
        for daypath in daylist:
            datelist+=[daypath.split('/')[len(daypath.split('/'))-2]+'/'+\
                        daypath.split('/')[len(daypath.split('/'))-1]]
        daylist=sorted(list(set(datelist))) # Get unique values
        if len(daylist) > 1:
            print 'You have collected data over multiple days - slacker'
            print 'Will run over: '+str(len(daylist))+' unique days'
            #sys.exit()
        if 'prevdaypath' in locals():
            del prevdaypath # Explicitly remove the previous daypath from locals
        daycount=0
        for daypath in daylist:
            yeardir=daypath.split('/')[len(daypath.split('/'))-2]   # Will be of the form Y2014
            daydir=daypath.split('/')[len(daypath.split('/'))-1]    # Will be of the form R201.01
            if daycount < len(daylist)-1:
                nextdaypath=daylist[daycount+1]
            print '\n'+daydir
            if defaults.rawconv=='True':
                st=obsread(defaults.outdir+'/*/'+yeardir+'/'+daydir+'/*.m')
            else:
                st=obsread(defaults.outdir+'/'+yeardir+'/'+daydir+'/*.m')
            print 'Merging data'
            try:
                st = st.detrend('simple')    # Detrend data before filling
                st.merge(fill_value=0)  # merge data, filling missing data with zeros -
                                        # allows for writing to multiplexed miniseed
            except:
                print 'Could not merge data for this day - same IDs but different sampling rates likely'
                samp_rate=st[0].stats.sampling_rate
                if 'st_dummy' in locals():
                    del st_dummy
                for tr in st:
                    if not tr.stats.sampling_rate==samp_rate:
                        print 'station: '+tr.stats.station+' samp-rate: '+\
                                str(tr.stats.sampling_rate)
                    else:
                        if 'st_dummy' in locals():
                            st_dummy+=tr
                        else:
                            st_dummy=Stream(tr)
                st=st_dummy
                st = st.detrend('simple')
                st.merge(fill_value=0)
            if defaults.debug==1:
                print 'I have read in '+str(len(st))+' traces'
                print 'They start at: '+str(st[0].stats.starttime.year)+'/'+\
                        str(st[0].stats.starttime.month)+'/'+\
                        str(st[0].stats.starttime.day)+' '+\
                        str(st[0].stats.starttime.hour)+':'+\
                        str(st[0].stats.starttime.minute)+':'+\
                        str(st[0].stats.starttime.second)
            if not defaults.alltime:     # If we are not running through all the
                                         # then there should be data for the
                                         # previous day to read in
                try:
                    prevday=st[0].stats.starttime-86400
                    prevdaypath=glob.glob(defaults.outdir+'/Y'+str(prevday.year)+
                                          '/R'+str(prevday.julday)+'.01')[0]
                except:
                    print 'Unable to find directory containing previous days data'
            # Work out true starttime and change channel names
            datastart=st[0].stats.starttime
            dataend=st[0].stats.endtime
            if daycount==0:
                startdate=datastart
            if daycount==len(daylist)-1:
                enddate=datastart
            for tr in st:
                if tr.stats.starttime < datastart:
                    datastart=tr.stats.starttime
                if tr.stats.endtime > dataend:
                    dataend=tr.stats.endtime
            # Read in more data once we have the start and end times well defined
            # Look for files in the day before and day after to get a complete file
            if 'prevdaypath' in locals():
                prevyeardir=prevdaypath.split('/')[len(prevdaypath.split('/'))-2]
                prevdaydir=prevdaypath.split('/')[len(prevdaypath.split('/'))-1]
                if defaults.rawconv:
                    if defaults.debug==1:
                        print 'Trying to read previous data from: '+defaults.outdir+'/*/'+\
                                prevyeardir+'/'+prevdaydir+'/*.*.23.*.m'
                    st+=obsread(defaults.outdir+'/*/'+prevyeardir+'/'+prevdaydir+'/*.*.23.*.m')
                else:
                    if defaults.debug==1:
                        print 'Trying to read previous data from: '+defaults.outdir+'/'+\
                                prevyeardir+'/'+prevdaydir+'/*.*.23.*.m'
                    st+=obsread(defaults.outdir+'/'+prevyeardir+'/'+prevdaydir+'/*.*.23.*.m')

                try:
                    st = st.detrend('simple')
                    st.merge(fill_value=0)  # merge data, filling missing data with zeros -
                                        # allows for writing to multiplexed miniseed
                except:
                    print 'Could not merge data for this day - same IDs but different sampling rates likely'
                    samp_rate=st[0].stats.sampling_rate
                    if 'st_dummy' in locals():
                        del st_dummy
                    for tr in st:
                        if not tr.stats.sampling_rate==samp_rate:
                            print 'station: '+tr.stats.station+' samp-rate: '+\
                                    str(tr.stats.sampling_rate)
                        else:
                            if 'st_dummy' in locals():
                                st_dummy+=tr
                            else:
                                st_dummy=Stream(tr)
                    st=st_dummy
                    st = st.detrend('simple')
                    st.merge(fill_value=0)



                if defaults.debug==1:
                    print 'I have read in '+str(len(st))+' traces'
                    print 'They start at: '+str(st[0].stats.starttime.year)+'/'+\
                            str(st[0].stats.starttime.month)+'/'+\
                            str(st[0].stats.starttime.day)+' '+\
                            str(st[0].stats.starttime.hour)+':'+\
                            str(st[0].stats.starttime.minute)+':'+\
                            str(st[0].stats.starttime.second)
            if 'nextdaypath' in locals():
                nextyeardir=nextdaypath.split('/')[len(nextdaypath.split('/'))-2]
                nextdaydir=nextdaypath.split('/')[len(nextdaypath.split('/'))-1]
                if defaults.rawconv:
                    if defaults.debug==1:
                        print 'Trying to read next data from: '+defaults.outdir+'/*/'+\
                                nextyeardir+'/'+nextdaydir+'/*.*.00.*.m'
                    st+=obsread(defaults.outdir+'/*/'+nextyeardir+'/'+nextdaydir+'/*.*.00.*.m')
                else:
                    if defaults.debug==1:
                        print 'Trying to read next data from: '+defaults.outdir+'/'+\
                                nextyeardir+'/'+nextdaydir+'/*.*.00.*.m'
                    st+=obsread(defaults.outdir+'/'+nextyeardir+'/'+nextdaydir+'/*.*.00.*.m')
                try:
                    st.detrend('simple')
                    st.merge(fill_value=0) # merge data filling gaps
                    #st.merge(fill_value=0)  # merge data, filling missing data with zeros -
                                            # allows for writing to multiplexed miniseed
                except:
                    print 'Could not merge data for this day - same IDs but different sampling rates likely'
                    samp_rate=st[0].stats.sampling_rate
                    if 'st_dummy' in locals():
                        del st_dummy
                    for tr in st:
                        if not tr.stats.sampling_rate==samp_rate:
                            print 'station: '+tr.stats.station+' samp-rate: '+\
                                    str(tr.stats.sampling_rate)
                        else:
                            if 'st_dummy' in locals():
                                st_dummy+=tr
                            else:
                                st_dummy=Stream(tr)
                    st=st_dummy
                    st.detrend('simple')
                    st.merge(fill_value=0)

                if defaults.debug==1:
                    print 'I have read in '+str(len(st))+' traces'
                    print 'They start at: '+str(st[0].stats.starttime.year)+'/'+\
                            str(st[0].stats.starttime.month)+'/'+\
                            str(st[0].stats.starttime.day)+' '+\
                            str(st[0].stats.starttime.hour)+':'+\
                            str(st[0].stats.starttime.minute)+':'+\
                            str(st[0].stats.starttime.second)

##########################Download GeoNet data for day if required############
            if defaults.getGeoNet:
                print 'Downloading GeoNet data for the day'
                oldlen = len(st)
                import subprocess
                pad=''
                for station in defaults.geostalist:
                    if len(station.name) == 3:
                        pad = '..'
                    elif len(station.name) == 4:
                        pad = '.'
                    for channel in station.channels:
                        subprocess.call(['java','-jar','pro/GeoNetCWBQuery-4.2.0-bin.jar',
                            '-d','1d','-t','ms','-b',day.year+day.julday+' 00:00:00',
                            '-s','NZ'+station.name+pad+channel,
                            '-o',daypath+'/%z%y%M%D.%s.%c.ms'])
                        st+=read(daypath+'/*.'+station.name+'.'+channel+'.ms')
                print 'I have downloaded and read in an extra '+\
                        str(len(st)-oldlen)+' traces'

            prevdaypath=daypath
            # Change the channel names from those output by rt2ms (101,102,103) to the
            # true channel names taken from the station definitions above
            if defaults.rename:
                for tr in st:
                    ksta=0
                    for sta in defaults.stalist:
                        if tr.stats.station == sta.name:
                            staid=ksta
                        else:
                            ksta+=1
                    if tr.stats.channel == '101':
                        tr.stats.channel=defaults.stalist[staid].channels[0]
                    if tr.stats.channel == '102':
                        tr.stats.channel=defaults.stalist[staid].channels[1]
                    if tr.stats.channel == '103':
                        tr.stats.channel=defaults.stalist[staid].channels[3]


            # Archive data here
            if defaults.archive:
                for tr in st:
                    daystart=UTCDateTime(yeardir[1:5]+daydir[1:4])
                    # print daystart
                    # print daystart+86400
                    tr=tr.trim(daystart,daystart+86400)
                    Ypath=defaults.arcdir+'/'+yeardir
                    Dpath=Ypath+'/'+daydir
                    if not os.path.isdir(Ypath):
                        os.mkdir(Ypath)
                    if not os.path.isdir(Dpath):
                        os.mkdir(Dpath)
                    outpath=Dpath+'/'+tr.stats.station+'.'+tr.stats.network+'..'+\
                            tr.stats.channel+'.'+str(tr.stats.starttime.year)+'.'+\
                            str(tr.stats.starttime.julday).zfill(3)
                    if len(tr.data) > 0:
                        tr.write(outpath,format="MSEED", encoding=defaults.dformat)
                    sys.stdout.write('Archived file written as: '+outpath+'\r')
                    sys.stdout.flush()
            print '\n'
            # If data are longer than your desired file length then cut
            if defaults.merge:
                datalen=dataend-datastart
                if datalen > defaults.filelenwant:
                    for hour in range(1,int(round(datalen/defaults.filelenwant)+1)):
                        st1=deepcopy(st)
                        # st1=st.copy() # Copy file to preserve the original data
                        cutstart=datastart+((hour-1)*defaults.filelenwant)
                        # Edit cut start to be cutting at the hour
                        cutstart=UTCDateTime(cutstart.year,cutstart.month,\
                                             cutstart.day,cutstart.hour)
                        cutend=cutstart+defaults.filelenwant
                        st1=st1.trim(cutstart,cutend)
                        filename=defaults.contbase+'/'+str(st1[0].stats.starttime.year)+'/'+\
                                str(st1[0].stats.starttime.month).zfill(2)+'/'+\
                                str(st1[0].stats.starttime.year)+'-'+\
                                str(st1[0].stats.starttime.month).zfill(2)+'-'+\
                                str(st1[0].stats.starttime.day).zfill(2)+'-'+\
                                str(st1[0].stats.starttime.hour).zfill(2)+\
                                str(st1[0].stats.starttime.minute).zfill(2)+'-'+\
                                str(st1[0].stats.starttime.second).zfill(2)+'.'+\
                                defaults.network+'_'+str(len(st1)).zfill(3)+'_00'
                        # Now write
                        try:
                            if not os.path.isdir(defaults.contbase+'/'+\
                                            str(st1[0].stats.starttime.year)):
                                os.mkdir(defaults.contbase+'/'+\
                                           str(st1[0].stats.starttime.year))
                            if not os.path.isdir(defaults.contbase+'/'+\
                                            str(st1[0].stats.starttime.year)+'/'+\
                                            str(st1[0].stats.starttime.month).zfill(2)):
                                os.mkdir(defaults.contbase+'/'+\
                                            str(st1[0].stats.starttime.year)+'/'+\
                                            str(st1[0].stats.starttime.month).zfill(2))
                            st1.write(filename,format="MSEED",endcoding=defaults.dformat)
                            sys.stdout.write('File written as: '+filename+'\r')
                            sys.stdout.flush()
                        except:
                            print('File not written, unknown error, check log:')
                            print(filename)
                            for tr in st1:
                                f.write(filename+'\n')
                                f.write(str(tr.stats))
                                f.write('\n\n')
                else:
                    filename=defaults.contbase+'/'+str(st[0].stats.starttime.year)+'/'+\
                        str(st[0].stats.starttime.month).zfill(2)+'/'+\
                        str(st[0].stats.starttime.year)+'-'+\
                        str(st[0].stats.starttime.month).zfill(2)+'-'+\
                        str(st[0].stats.starttime.day).zfill(2)+'-'+\
                        str(st[0].stats.starttime.hour).zfill(2)+\
                        str(st[0].stats.starttime.minute).zfill(2)+'-'+\
                        str(st[0].stats.starttime.second).zfill(2)+'.'+\
                        defaults.network+'_'+str(len(st)).zfill(3)+'_00'
                    st.write(filename,format="MSEED",endcoding=defaults.dformat)
                    sys.stdout.write('File written as: '+filename+'\r')
                    sys.stdout.flush()
            st=[]
            daycount+=1

        # END of mergin section, you should now have a continuous database!

        # Remove non-merged data to conserve space
        if defaults.rmold:
            rmlist=glob.glob(defaults.outdir+'/*/Y*')
            for tree in rmlist:
                print 'Removing old data'
                shutil.rmtree(tree)
    elif not defaults.detrout=='none':
        # Set up variables required for later
        print "\nAs I'm not merging anything, can you tell me the start date for the"
        print "data to detect over in the form yyyy/mm/dd?"
        startdate = raw_input('start date: > ')
        # Check input
        try:
            startdate=datetime.datetime.strptime(startdate,'%Y/%m/%d')
        except:
            print 'startdate not accepted, please check formatting'
            sys.exit()
        print "Can you now tell me the end date, in the same format?"
        enddate = raw_input('end date: > ')
        # Check input
        try:
            enddate=datetime.datetime.strptime(enddate,'%Y/%m/%d')
        except:
            print 'enddate not accepted, please check formatting'
            sys.exit()
        print "Startdate: "+str(startdate)
        print "Enddate: "+str(enddate)

    ################## Detection routine and clipping to triggered database########

    # Run condet
    if defaults.detrout=='condet':
        print('I will register and detect a month at a time - you will need to')
        print('respond to prompts')
        if startdate.year!=enddate.year:
            for regyear in range(startdate.year,enddate.year):
                for regmonth in range(startdate.month):
                    print 'bob'
                    call(["dirf",defaults.contbase+'/'+regyear+'/'+regmonth+'/*'])
                    call(["autoreg"])
                    call(["condet"])
    elif defaults.detrout=='obspy':
        # Run Calum's detection routine derived from obspy
        import pro.STALTA as STALTA
        startdate=str(startdate.year)+'/'+str(startdate.month).zfill(2)+'/'+\
                str(startdate.day).zfill(2)
        enddate=str(enddate.year)+'/'+str(enddate.month).zfill(2)+'/'+\
                str(enddate.day).zfill(2)
        if defaults.debug==1:
            print "STALTA.cjc_trigger_routine("+startdate+','+enddate+','+defaults.contbase\
                    +','+defaults.trigbase+','+defaults.routype+')'
        wavelist=STALTA.cjc_trigger_routine(startdate,enddate,defaults.contbase,\
                defaults.trigbase,defaults.routype)
        # Generate file list which can then be autoregistered
        # f1.open('filenr.lis','w')
        # listno=0
        # for wavefile in wavelist:
            # listno+=1
            # f1.write('#'+str(listno).rjust(3)+'  '+wavefile+'\n')
        # f1.close()
        print 'Generating s-files'
        import pro.Sfile_util as Sfile_util
        import ntpath
        sfilelist=[]
        for wavepath in wavelist:
            wavefile=ntpath.basename(wavepath)
            shutil.copy(wavepath,wavefile)   # Copy the file locally so that we can
                                            # cope with filename concatenation
            sfilelist.append(Sfile_util.blanksfile(wavefile,'L',defaults.userID,'.',defaults.overwrite))
            # Write to local directory, otherwise fortran concatenation errors are
            # likely in the filename when running the filterdefaults.picker routine - move
            # later in the script
    else:
        # wavelist=glob.glob('*-*-*-*-*.*_*_*')
        sfilelist=glob.glob('*L.S*')
        import ntpath
        print('No detection routine selected, you just have continuous data')

    ###################Picking routines############################################
    if defaults.picker=='FP':
        print('Will now run the filterpicker routine')
        # run filter defaults.picker routine adapted within rtquake and altered for use here
        # by Calum Chamberlain.
        i=0
        for sfile in sfilelist:
            wavefile=Sfile_util.readwavename(sfile)
            # wavefile=ntpath.basename(wavelist[i])
            if not os.path.isfile(wavefile):
                print 'Wavefile '+wavefile+' not found, will not locate'
                continue
            i+=1
            print(["pro/rdtrigL","-sfile",sfile,"-wavefile",wavefile,\
                    "-locate","1","-maxres",str(defaults.maxres),\
                    "-keep","1"])
            call(["pro/rdtrigL","-sfile",sfile,"-wavefile",wavefile,\
                    "-locate","1","-maxres",str(defaults.maxres),\
                    "-keep","1"])
            if not os.path.isdir(defaults.sfilebase+'/'+sfile[13:17]):
                os.mkdir(defaults.sfilebase+'/'+sfile[13:17])
            if not os.path.isdir(defaults.sfilebase+'/'+sfile[13:17]+'/'\
                                 +sfile[17:19]):
                os.mkdir(defaults.sfilebase+'/'+sfile[13:17]+'/'\
                    +sfile[17:19])
            shutil.move(sfile,defaults.sfilebase+'/'+sfile[13:17]+'/'\
                    +sfile[17:19]+'/'+sfile) # Move the picked sfile to the
                                                   # correct place.
            print 'Written file as: '+sfile,defaults.sfilebase+'/'+sfile[13:17]+'/'\
                    +sfile[17:19]+'/'+sfile
            os.remove(wavefile) # Remove local wavefiles
    f.close()


