#!/usr/bin/python

# Code to generate an empty s-file in nordic format for a given input waveform
# file.  Use to replace seisan's autoreg with a tool that doesn't require
# interactive responses

# v0.1 - Beginning writing in Paris-Orly airport while very sleepy
# v1.0 - Finished and working, impliments everything desired at this point

# Calum John Chamberlain, Victoria University of Wellington - 11/11/2014

def blanksfile(wavefile,evtype,userID,outdir,overwrite):

###############################################################################

    # Arguments are the path of a wavefile (multiplexed miniseed file required)
    # Event type (L,R,D) and user ID (four characters as used in seisan)

###############################################################################

    # Example s-file format:
    # 2014  719  617 50.2 R                                                         1
    # ACTION:ARG 14-11-11 10:53 OP:CALU STATUS:               ID:20140719061750     I
    # 2014/07/2014-07-19-0617-50.SAMBA_030_00                                       6
    # STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7

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
