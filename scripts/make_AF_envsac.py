#!/usr/bin/python
"""
Simple function to read in day long miniseed data and generate 1Hz bandpassed
envelopes for the use in tremor detection.  Based on matlab code of Aaron Wech.
"""

def envsac(tr, low, high, delta, debug=False):
    """
    Function to generate an envelope from given seismic data - requires obspy

    :type tr: obspy.Trace
    :type low: float
    :param low: Lowcut in Hz
    :type high: float
    :param high: Highcut in Hz
    :type delta: float
    :param delta: Sampling interval desired

    :return env: Envelope as an obspy.Trace object
    """
    from obspy.signal.filter import bandpass, envelope
    if debug:
        print 'Detrend'
    env=tr.detrend('simple')                    # Demean the data
    del tr
    if debug:
        print 'Resample'
    env.resample(40)                            # Downsample to 40Hz
    if debug:
        print 'Bandpass'
    env.data=bandpass(env.data, low, high, env.stats.sampling_rate,\
                      3, True)                  # Filter
    if debug:
        print 'Envelope'
    env.data=envelope(env.data)                 # Generate envelope
    if debug:
        print 'Resample'
    env.resample(1)                             # Deimate to 1 Hz
    return env

def check_daylong(tr):
    """
    Function to check the data quality of the daylong file - check to see that
    the day isn't just zeros, with large steps, if it is then the resampling will
    hate it.

    :type tr: obspy.Trace

    :return qual: bool
    """
    import numpy as np
    if len(tr.data)-len(np.nonzero(tr.data)) < 0.5*len(tr.data):
        qual=False
    else:
        qual=True
    return qual

if __name__ == '__main__':
    debug=True
    debugplot=False
    import sys, os
    if not len(sys.argv) == 5:
        raise IOError('Insufficent arguments, I need startdate, enddate, indir and outdir')
    else:
        startdate=str(sys.argv[1])
        enddate=str(sys.argv[2])
        indir=str(sys.argv[3])
        outdir=str(sys.argv[4])
    from obspy import UTCDateTime
    startdate=UTCDateTime(startdate[0:4]+'-'+startdate[4:6]+'-'+startdate[6:8])
    enddate=UTCDateTime(enddate[0:4]+'-'+enddate[4:6]+'-'+enddate[6:8])
    kdays=int((enddate-startdate)/86400)+1
    print 'I will loop through '+str(kdays)+' days'
    import glob
    from obspy import read as obsread
    for i in xrange(kdays):
        date=startdate+86400*i
        print 'Working on day: '+str(date)
        daydir=str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)
        infiles=glob.glob(indir+'/Y'+str(date.year)+'/R'+str(date.julday).zfill(3)+\
                         '.01/*')
        for infile in infiles:
            if debug:
                print 'Reading in '+infile
            tr=obsread(infile)
            # Fill any gaps in the data
            tr=tr.merge(fill_value='interpolate')
            # Make daylong
            tr=tr.detrend('simple')
            tr=tr.trim(starttime=date, endtime=date+86400, pad=True, fill_value=0,\
                    nearest_sample=False)
            tr=tr[0]
            if debug:
                print 'Read in file, it is '+str(len(tr.data))+' samples long'
            qual=check_daylong(tr)
            if qual:
                trenv=envsac(tr, 1.5, 5.0, 1, debug)
                if debugplot:
                    trenv.plot()
                del tr
                if not os.path.isdir(outdir+'/'+daydir):
                    os.makedirs(outdir+'/'+daydir)
                trenv.write(outdir+'/'+daydir+'/'+daydir+'.'+trenv.stats.network+'.'+\
                            trenv.stats.station+'.'+trenv.stats.channel+'_env.sac',\
                            format='SAC')
                print 'Written file as: '+outdir+'/'+daydir+'/'+daydir+'.'+\
                        trenv.stats.network+'.'+trenv.stats.station+'.'+\
                        trenv.stats.channel+'_env.sac'
                del trenv
            else:
                msg="Trace contains more zeros than data, check this!" + \
                        infile
                raise ValueError(msg)
                del tr
