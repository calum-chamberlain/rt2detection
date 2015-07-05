#!/usr/bin/python
"""
Part of the rt2detection module to read nordic format s-files and write them
EQcorrscan is a python module designed to run match filter routines for
seismology, within it are routines for integration to seisan and obspy.
With obspy integration (which is necessary) all main waveform formats can be
read in and output.

Code generated by Calum John Chamberlain of Victoria University of Wellington,
2015.

All rights reserved.

"""

# Import PICK class definition from makeSfile

from obspy import UTCDateTime
class PICK:
    """
    Pick information for seisan implimentation
    """
    pickcount=0
    def __init__(self, station, channel, impulsivity, phase, weight, polarity,
                 time, coda, amplitude, peri, azimuth, velocity, AIN, SNR,
                 azimuthres, timeres, finalweight, distance, CAZ):
        self.station=station
        self.channel=channel
        self.impulsivity=impulsivity
        self.phase=phase
        self.weight=weight
        self.polarity=polarity
        self.time=time
        self.coda=coda
        self.amplitude=amplitude
        self.peri=peri
        self.azimuth=azimuth
        self.velocity=velocity
        self.AIN=AIN
        self.SNR=SNR
        self.azimuthres=azimuthres
        self.timeres=timeres
        self.finalweight=finalweight
        self.distance=distance
        self.CAZ=CAZ
        self.pickcount+=1

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

def _int_conv(string):
    """
    Convenience tool to convert from string to integer, if empty string return
    a 999 rather than an error
    """
    try:
        intstring=int(string)
    except:
        intstring=999
    return intstring

def _float_conv(string):
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
            if int(topline[16:18]) == 60:
                sfilename_header.time=UTCDateTime(int(topline[1:5]),int(topline[6:8]),
                                        int(topline[8:10]),int(topline[11:13]),
                                        int(topline[13:15])+1,0
                                        ,int(topline[19:20])*10)
            else:
                sfilename_header.time=UTCDateTime(int(topline[1:5]),int(topline[6:8]),
                                        int(topline[8:10]),int(topline[11:13]),
                                        int(topline[13:15]),int(topline[16:18])
                                        ,int(topline[19:20])*10)
        except:
            sfilename_header.time=UTCDateTime(0)
        sfilename_header.loc_mod_ind=topline[20]
        sfilename_header.dist_ind=topline[21]
        sfilename_header.ev_id=topline[22]
        sfilename_header.latitude=_float_conv(topline[23:30])
        sfilename_header.longitude=_float_conv(topline[31:38])
        sfilename_header.depth=_float_conv(topline[39:43])
        sfilename_header.depth_ind=topline[44]
        sfilename_header.loc_ind=topline[45]
        sfilename_header.agency=topline[46:48].strip()
        sfilename_header.nsta=_int_conv(topline[49:51])
        sfilename_header.t_RMS=_float_conv(topline[52:55])
        sfilename_header.Mag_1=_float_conv(topline[56:59])
        sfilename_header.Mag_1_type=topline[60]
        sfilename_header.Mag_1_agency=topline[61:63].strip()
        sfilename_header.Mag_2=_float_conv(topline[64:67])
        sfilename_header.Mag_2_type=topline[68]
        sfilename_header.Mag_2_agency=topline[69:71].strip()
        sfilename_header.Mag_3=_float_conv(topline[72:75])
        sfilename_header.Mag_3_type=topline[76].strip()
        sfilename_header.Mag_3_agency=topline[77:79].strip()
    else:
        for line in f:
            if line[79]=='1':
                line=topline
                try:
                    sfilename_header.time=UTCDateTime(int(topline[1:5]),int(topline[6:8]),
                                                int(topline[8:10]),int(topline[11:13]),
                                                int(topline[13:15]),int(topline[16:18])
                                                ,int(topline[19:20])*10)
                except:
                    sfilename_header.time=UTCDateTime(0)
                sfilename_header.loc_mod_ind=topline[21]
                sfilename_header.dist_ind=topline[22]
                sfilename_header.ev_id=topline[23]
                sfilename_header.latitude=_float_conv(topline[24:30])
                sfilename_header.longitude=_float_conv(topline[31:38])
                sfilename_header.depth=_float_conv(topline[39:43])
                sfilename_header.depth_ind=topline[44]
                sfilename_header.loc_ind=topline[45]
                sfilename_header.agency=topline[46:48].strip()
                sfilename_header.nsta=_int_conv(topline[49:51])
                sfilename_header.t_RMS=_float_conv(topline[52:55])
                sfilename_header.Mag_1=_float_conv(topline[56:59])
                sfilename_header.Mag_1_type=topline[60]
                sfilename_header.Mag_1_agency=topline[61:63].strip()
                sfilename_header.Mag_2=_float_conv(topline[64:67])
                sfilename_header.Mag_2_type=topline[68]
                sfilename_header.Mag_2_agency=topline[69:71].strip()
                sfilename_header.Mag_3=_float_conv(topline[72:75])
                sfilename_header.Mag_3_type=topline[76].strip()
                sfilename_header.Mag_3_agency=topline[77:79].strip()
            if line[79]=='7':
                break
    f.close()
    return sfilename_header

def readpicks(sfilename):
    """
    Function to read pick informaiton from the s-file

    :type sfilename: String

    :return: Sfile_tile.PICK
    """
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
        try:
            time=UTCDateTime(evtime.year,evtime.month,evtime.day,
                             int(line[18:20]),int(line[20:22]),int(line[23:25]),
                             int(line[26:28]))
        except (ValueError):
            time=UTCDateTime(evtime.year,evtime.month,evtime.day,
                             int(line[18:20]),int(line[20:22]),0,0)
            time+=60 # Add 60 seconds on to the time, this copes with s-file
        coda=_int_conv(line[28:33])
        amplitude=_float_conv(line[34:40])
        peri=_float_conv(line[41:45])
        azimuth=_float_conv(line[46:51])
        velocity=_float_conv(line[52:56])
        if header[57:60]=='AIN':
            SNR=''
            AIN=_int_conv(line[57:60])
        elif header[57:60]=='SNR':
            AIN=''
            SNR=_float_conv(line[57:60])
        azimuthres=_int_conv(line[60:63])
        timeres=_float_conv(line[63:70])
        finalweight=_int_conv(line[70])
        distance=_float_conv(line[71:75])
        CAZ=_int_conv(line[76:79])
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
    wavename=[]
    for line in f:
        if line[79]=='6':
            wavename.append(line[1:79].strip())
    f.close()
    return wavename


def blanksfile(wavefile,evtype,userID,outdir,overwrite):
    """
    Module to generate an empty s-file with a populated header for a given
    waveform.

###############################################################################

    # Arguments are the path of a wavefile (multiplexed miniseed file required)
    # Event type (L,R,D) and user ID (four characters as used in seisan)

###############################################################################

    # Example s-file format:
    # 2014  719  617 50.2 R                                                         1
    # ACTION:ARG 14-11-11 10:53 OP:CALU STATUS:               ID:20140719061750     I
    # 2014/07/2014-07-19-0617-50.SAMBA_030_00                                       6
    # STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7
    """

    from obspy import read as obsread
    import sys,os, datetime
    try:
        st=obsread(wavefile)
    except:
        print 'Wavefile: '+wavefile+' is invalid, try again with real data.'
        sys.exit()
    # Check that user ID is the correct length
    if len(userID) != 4:
        print 'User ID must be 4 characters long'
        sys.exit()
    # Check that outdir exists
    if not os.path.isdir(outdir):
        print 'Out path does not exist, I will not create this: '+outdir
        sys.exit()
    # Check that evtype is one of L,R,D
    if evtype not in ['L','R','D']:
        print 'Event type must be either L, R or D'
        sys.exit()

    # Generate s-file name in the format dd-hhmm-ss[L,R,D].Syyyymm
    sfilename=outdir+'/'+str(st[0].stats.starttime.day).zfill(2)+'-'+\
            str(st[0].stats.starttime.hour).zfill(2)+\
            str(st[0].stats.starttime.minute).zfill(2)+'-'+\
            str(st[0].stats.starttime.second).zfill(2)+evtype+'.S'+\
            str(st[0].stats.starttime.year)+\
            str(st[0].stats.starttime.month).zfill(2)
    # Check is sfilename exists
    if os.path.isfile(sfilename) and overwrite=='False':
        print 'Desired sfilename: '+sfilename+' exists, will not overwrite'
        for i in range(1,10):
            sfilename=outdir+'/'+str(st[0].stats.starttime.day).zfill(2)+'-'+\
                    str(st[0].stats.starttime.hour).zfill(2)+\
                    str(st[0].stats.starttime.minute).zfill(2)+'-'+\
                    str(st[0].stats.starttime.second+i).zfill(2)+evtype+'.S'+\
                    str(st[0].stats.starttime.year)+\
                    str(st[0].stats.starttime.month).zfill(2)
            if not os.path.isfile(sfilename):
                break
        else:
            print 'Tried generated files up to 10s in advance and found they'
            print 'all exist, you need to clean your stuff up!'
            sys.exit()
        # sys.exit()
    f=open(sfilename,'w')
    # Write line 1 of s-file
    f.write(' '+str(st[0].stats.starttime.year)+' '+\
            str(st[0].stats.starttime.month).rjust(2)+\
            str(st[0].stats.starttime.day).rjust(2)+' '+\
            str(st[0].stats.starttime.hour).rjust(2)+\
            str(st[0].stats.starttime.minute).rjust(2)+' '+\
            str(st[0].stats.starttime.second).rjust(4)+' '+\
            evtype+'1'.rjust(58)+'\n')
    # Write line 2 of s-file
    f.write(' ACTION:ARG '+str(datetime.datetime.now().year)[2:4]+'-'+\
            str(datetime.datetime.now().month).zfill(2)+'-'+\
            str(datetime.datetime.now().day).zfill(2)+' '+\
            str(datetime.datetime.now().hour).zfill(2)+':'+\
            str(datetime.datetime.now().minute).zfill(2)+' OP:'+\
            userID.ljust(4)+' STATUS:'+'ID:'.rjust(18)+\
            str(st[0].stats.starttime.year)+\
            str(st[0].stats.starttime.month).zfill(2)+\
            str(st[0].stats.starttime.day).zfill(2)+\
            str(st[0].stats.starttime.hour).zfill(2)+\
            str(st[0].stats.starttime.minute).zfill(2)+\
            str(st[0].stats.starttime.second).zfill(2)+\
            'I'.rjust(6)+'\n')
    # Write line 3 of s-file
    f.write(' '+wavefile+'6'.rjust(79-len(wavefile))+'\n')
    # Write final line of s-file
    f.write(' STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU'+\
            ' VELO AIN AR TRES W  DIS CAZ7\n')
    f.close()
    print 'Written s-file: '+sfilename
    return sfilename

def populateSfile(sfilename, picks):
    """
    Module to populate a blank nordic format S-file with pick information,
    arguments required are the filename of the blank s-file and the picks
    where picks is a dictionary of picks including station, channel,
    impulsivity, phase, weight, polarity, time, coda, amplitude, peri, azimuth,
    velocity, SNR, azimuth residual, Time-residual, final weight,
    epicentral distance & azimuth from event to station.

    This is a full pick line information from the seisan manual, P. 341
    """

    f=open(sfilename, 'r')
    # Find type 7 line, under which picks should be - if there are already
    # picks there we should preserve them
    lineno=0
    body=''
    header=''
    if 'headerend' in locals():
        del headerend
    for line in f:
        identifier=line[79]
        if 'headerend' in locals():
            body+=line
        else:
            header+=line
        if identifier=='7':
            headerend=lineno
        lineno+=1
    f.close()
    #
    # Now generate lines for the new picks
    newpicks=''
    for pick in picks:
        newpicks+=' '+pick.station.ljust(5)+pick.channel[0]+\
                pick.channel[len(pick.channel)-1]+' '+pick.impulsivity+\
                pick.phase.ljust(4)+str(pick.weight).rjust(1)+' '+\
                pick.polarity+' '+str(pick.time.hour).rjust(2)+\
                str(pick.time.minute).rjust(2)+str(pick.time.second).rjust(3)+\
                '.'+str(pick.time.microsecond).ljust(2)+\
                str(pick.coda).rjust(5)+str(pick.amplitude).rjust(7)+\
                str(pick.peri).rjust(5)+str(pick.azimuth).rjust(6)+\
                str(pick.velocity).rjust(5)+str(pick.AIN).rjust(4)+\
                str(pick.azimuthres).rjust(3)+str(pick.timeres).rjust(5)+\
                str(pick.finalweight).rjust(2)+str(pick.distance).rjust(4)+\
                str(pick.CAZ).rjust(4)+' \n'
    # Write all new and old info back in
    f=open(sfilename, 'w')
    f.write(header)
    f.write(newpicks)
    f.write(body)
    f.close()

if __name__=='__main__':
    # Read arguments
    import sys, os
    if len(sys.argv) != 6:
        print 'Requires 5 arguments: wavefile, evtype, userID, outdir, overwrite'
        sys.exit()
    else:
        wavefile=str(sys.argv[1])
        evtype=str(sys.argv[2])
        userID=str(sys.argv[3])
        outdir=str(sys.argv[4])
        overwrite=str(sys.argv[5])
    sfilename=blanksfile(wavefile,evtype,userID,outdir,overwrite)
    print sfilename
