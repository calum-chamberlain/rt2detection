#!/usr/bin/python
"""
Part of the EQcorrscan module to read nordic format s-files
EQcorrscan is a python module designed to run match filter routines for
seismology, within it are routines for integration to seisan and obspy.
With obspy integration (which is necessary) all main waveform formats can be
read in and output.

This main section contains a script, LFE_search.py which demonstrates the usage
of the built in functions from template generation from picked waveforms
through detection by match filter of continuous data to the generation of lag
times to be used for relative locations.

The match-filter routine described here was used a previous Matlab code for the
Chamberlain et al. 2014 G-cubed publication.  The basis for the lag-time
generation section is outlined in Hardebeck & Shelly 2011, GRL.

Code generated by Calum John Chamberlain of Victoria University of Wellington,
2015.

All rights reserved.

Pre-requisites:
    gcc             - for the installation of the openCV correlation routine
    python-joblib   - used for parallel processing
    python-obspy    - used for lots of common seismological processing
                    - requires:
                        numpy
                        scipy
                        matplotlib
    python-pylab    - used for plotting
"""

# Import PICK class definition from makeSfile
from obspy import UTCDateTime

class EVENTINFO:
    def __init__(self, time, loc_mod_ind, dist_ind, ev_id, latitude, longitude,
                 depth, depth_ind, loc_ind, agency, nsta, t_RMS, Mag_1,
                 Mag_1_type, Mag_1_agency, Mag_2, Mag_2_type, Mag_2_agency,
                 Mag_3, Mag_3_type, Mag_3_agency):
        self.time=time
        self.loc_mod_ind=loc_mod_ind
        self.dist_ind=dist_ind
        self.ev_id=ev_id
        self.latitude=latitude
        self.longitude=longitude
        self.depth=depth
        self.depth_ind=depth_ind
        self.loc_ind=loc_ind
        self.agency=agency
        self.nsta=nsta
        self.t_RMS=t_RMS
        self.Mag_1=Mag_1
        self.Mag_1_type=Mag_1_type
        self.Mag_1_agency=Mag_1_agency
        self.Mag_2=Mag_2
        self.Mag_2_type=Mag_2_type
        self.Mag_2_agency=Mag_2_agency
        self.Mag_3=Mag_3
        self.Mag_3_type=Mag_3_type
        self.Mag_3_agency=Mag_3_agency

def int_conv(string):
    """
    Convenience tool to convert from string to integer, if empty string return
    a 999 rather than an error
    """
    try:
        intstring=int(string)
    except:
        intstring=999
    return intstring

def float_conv(string):
    """
    Convenience tool to convert from string to float, if empty string return
    NaN rather than an error
    """
    try:
        floatstring=float(string)
    except:
        floatstring=float('NaN')
    return floatstring

def readheader(sfilename):
    f=open(sfilename,'r')
    sfilename_header=EVENTINFO(UTCDateTime(), '', '', '',  float('NaN'),
                               float('NaN'), float('NaN'), '', '', '', 0,
                               float('NaN'), float('NaN'), '', '', float('NaN'),
                               '', '', float('NaN'), '', '')
    topline=f.readline()
    if topline[79]==' ' or topline[79]=='1':
        # Topline contains event information
        try:
            sfilename_header.time=UTCDateTime(int(topline[1:5]),int(topline[6:8]),
                                        int(topline[9:11]),int(topline[11:13]),
                                        int(topline[13:15]),int(topline[16:18])
                                        ,int(topline[19:20])*10)
        except:
            sfilename_header.time=UTCDateTime(0)
        sfilename_header.loc_mod_ind=topline[21]
        sfilename_header.dist_ind=topline[22]
        sfilename_header.ev_id=topline[23]
        sfilename_header.latitude=float_conv(topline[24:30])
        sfilename_header.longitude=float_conv(topline[31:38])
        sfilename_header.depth=float_conv(topline[39:43])
        sfilename_header.depth_ind=topline[44]
        sfilename_header.loc_ind=topline[45]
        sfilename_header.agency=topline[46:48].strip()
        sfilename_header.nsta=int_conv(topline[49:51])
        sfilename_header.t_RMS=float_conv(topline[52:55])
        sfilename_header.Mag_1=float_conv(topline[56:59])
        sfilename_header.Mag_1_type=topline[60]
        sfilename_header.Mag_1_agency=topline[61:63].strip()
        sfilename_header.Mag_2=float_conv(topline[64:67])
        sfilename_header.Mag_2_type=topline[68]
        sfilename_header.Mag_2_agency=topline[69:71].strip()
        sfilename_header.Mag_3=float_conv(topline[72:75])
        sfilename_header.Mag_3_type=topline[76].strip()
        sfilename_header.Mag_3_agency=topline[77:79].strip()
    else:
        for line in f:
            if line[79]=='1':
                line=topline
                try:
                    sfilename_header.time=UTCDateTime(int(topline[1:5]),int(topline[6:8]),
                                                int(topline[9:11]),int(topline[11:13]),
                                                int(topline[13:15]),int(topline[16:18])
                                                ,int(topline[19:20])*10)
                except:
                    sfilename_header.time=UTCDateTime(0)
                sfilename_header.loc_mod_ind=topline[21]
                sfilename_header.dist_ind=topline[22]
                sfilename_header.ev_id=topline[23]
                sfilename_header.latitude=float_conv(topline[24:30])
                sfilename_header.longitude=float_conv(topline[31:38])
                sfilename_header.depth=float_conv(topline[39:43])
                sfilename_header.depth_ind=topline[44]
                sfilename_header.loc_ind=topline[45]
                sfilename_header.agency=topline[46:48].strip()
                sfilename_header.nsta=int_conv(topline[49:51])
                sfilename_header.t_RMS=float_conv(topline[52:55])
                sfilename_header.Mag_1=float_conv(topline[56:59])
                sfilename_header.Mag_1_type=topline[60]
                sfilename_header.Mag_1_agency=topline[61:63].strip()
                sfilename_header.Mag_2=float_conv(topline[64:67])
                sfilename_header.Mag_2_type=topline[68]
                sfilename_header.Mag_2_agency=topline[69:71].strip()
                sfilename_header.Mag_3=float_conv(topline[72:75])
                sfilename_header.Mag_3_type=topline[76].strip()
                sfilename_header.Mag_3_agency=topline[77:79].strip()
            if line[79]=='7':
                break
    f.close()
    return sfilename_header

def readpicks(sfilename):
    from makeSfile import PICK
    # First we need to read the header to get the timing info
    sfilename_header=readheader(sfilename)
    evtime=sfilename_header.time
    f=open(sfilename,'r')
    pickline=[]
    lineno=0
    if 'headerend' in locals():
        del headerend
    for line in f:
        if 'headerend' in locals():
            if line[79]==' ':
                pickline+=[line]
        if line[79]=='7':
            header=line
            headerend=lineno
        lineno+=1
    picks=[]
    for line in pickline:
        if line[18:28].strip()=='': # If line is empty miss it
            continue
        station=line[1:6].strip()
        channel=line[6:8].strip()
        impulsivity=line[9]
        weight=line[15]
        if weight=='_':
            phase=line[10:17]
            weight=''
            polarity=''
        else:
            phase=line[10:14].strip()
            polarity=line[6]
        time=UTCDateTime(evtime.year,evtime.month,evtime.day,
                         int(line[18:20]),int(line[20:22]),int(line[23:25]),
                         int(line[26:28]))
        coda=int_conv(line[28:33])
        amplitude=float_conv(line[34:40])
        peri=float_conv(line[41:45])
        azimuth=float_conv(line[46:51])
        velocity=float_conv(line[52:56])
        if header[57:60]=='AIN':
            SNR=''
            AIN=int_conv(line[57:60])
        elif header[57:60]=='SNR':
            AIN=''
            SNR=float_conv(line[57:60])
        azimuthres=int_conv(line[60:63])
        timeres=float_conv(line[63:70])
        finalweight=int_conv(line[70])
        distance=float_conv(line[71:75])
        CAZ=int_conv(line[76:79])
        picks+=[PICK(station, channel, impulsivity, phase, weight, polarity,
                 time, coda, amplitude, peri, azimuth, velocity, AIN, SNR,
                 azimuthres, timeres, finalweight, distance, CAZ)]
    f.close()
    return picks

def readwavename(sfilename):
    """
    Convenience function to extract the waveform filename from the s-file,
    returns a list of waveform names found in the s-file as multiples can
    be present.
    """
    f=open(sfilename)
    for line in f:
        if line[79]=='6':
            if 'wavename' in locals():
                wavename+=line[1:79].strip()
            else:
                wavename=line[1:79].strip()
    f.close()
    return wavename
