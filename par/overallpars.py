#!/usr/bin/python
"""
# A python style declaration of variables used in rt2detection.py
"""
############################Code definitions###################################
indir='/Volumes/Taranaki_01/data/boseca/SAMBA_mar09/' # Location of raw reftek data
outdir='/media/Elements/SAMBA_archive/SAMBA_mseed'# Don't change temporary directory
<<<<<<< HEAD
contbase="/home/sw/seisan/WAV/SAMCO"              # Location for continuous
                                                  # waveforms
trigbase="/home/sw/seisan/WAV/SAMTR"              # Location for triggered
                                                  # waveforms
sfilebase="/home/sw/seisan/REA/SAMTR"             # Location for triggered
                                                  # s-files
network='SAMBA' # Needs to be five characters, e.g. SAMBA or RTAF_ - pad with
=======
contbase="/home/sw/seisan/WAV/COSAC"              # Location for continuous
                                                  # waveforms
trigbase="/home/sw/seisan/WAV/COSAW"              # Location for triggered
                                                  # waveforms
sfilebase="/home/sw/seisan/REA/COSAW"             # Location for triggered
                                                  # s-files
network='COSA_' # Needs to be five characters, e.g. SAMBA or RTAF_ - pad with
>>>>>>> 4eab03d65c72a9dc259ba578e5b556d1ab05752a
                # an underscore
archive=True    # Set to True to turn archiving on, will make daylong miniseed
                # files in the archive directory in Yyyyy/Rddd.01 type directories
arcdir='/Volumes/GeoPhysics_09/users-data/chambeca/SAMBA_archive/day_volumes_S'# Directory to archive data to
detrout='obspy' # Can either use the obspy triggering routine, or the seisan
                # condet routine (condet is a little less simple to use and is
                # not yet implimented).
routype='classic' # Required for calum's python based routine, see that file
                  # for options
#detrout='condet'
filelenwant=3600# Length to output continuous data in in seconds
dformat='STEIM2'# Format to output miniseed to
rmold=False     # Boolean - if you want to remove old data then set this to
                # True, it will remove the miniseed files in the MS_data
                # directory before it starts, then clean the non-merged data
                # at the end.
picker='FP'     # Picker method, can be set to False to not run one - calls
                # external routines
rawconv=True    # Set to true to look for new data in 'indir'
<<<<<<< HEAD
merge=False     # Set to true to merge data in MS_data into multiplexed
=======
merge=True      # Set to true to merge data in MS_data into multiplexed
>>>>>>> 4eab03d65c72a9dc259ba578e5b556d1ab05752a
                # miniseed files of length 'filelenwant'
maxres=1.5      # If picker is used, this is the maximum residual to allow for
                # location
userID='CALU'   # Your userID in seisan, must be four characters.
# userID='DFDP'
debug=1         # Set to 0 for a less verbose output
overwrite=False # Set to false to preserve old s-files, you will have to move
                # them manually in this case after the picking routine.
neo='False'     # If you have downloaded the data using neo then we will use
                # the routine for that. Currently set to string to allow for
                # a special 'Sometimes' case
<<<<<<< HEAD
rename=False    # Set to True to use the channel names defined below, if false,
=======
rename=True     # Set to True to use the channel names defined below, if false,
>>>>>>> 4eab03d65c72a9dc259ba578e5b556d1ab05752a
                # channels will be named as output by rt2ms (which at this stage
                # if called from here will be the wired 101,102,103 names.
                # I suggest only setting this to False if you have already
                # converted your raw data

# Descide whether you want it to run for all time or from a start to end date
<<<<<<< HEAD
alltime=False   # Set to true to run through all time, will desregard start and
                # end time arguments
startD='2014/03/26' # Start time for alltime=False in yyyy/mm/dd
endD = '2015/03/10' # End time for alltime=False in yyyy/mm/dd
getGeoNet =False# Boolean to get geonet data, if true requires the geostalist
=======
alltime=True    # Set to true to run through all time, will desregard start and
                # end time arguments
startD='2014/03/26' # Start time for alltime=False in yyyy/mm/dd
endD = '2015/03/10' # End time for alltime=False in yyyy/mm/dd
getGeoNet =True # Boolean to get geonet data, if true requires the geostalist
>>>>>>> 4eab03d65c72a9dc259ba578e5b556d1ab05752a
                # variable to be complete.
# Define local class
class STATION:
    'Station information for seismic station'
    stacount=0
    def __init__(self,name, netcode, loccode, das, channels):
        self.name=name
        self.netcode=netcode
        self.loccode=loccode
        self.das=das
        self.channels=channels
        self.stacount+=1

#########################USER EDIT SECTION-SETUP YOUR DATABASE#################

# Build station database: name, netcode, loccode, das, 101, 102, 103 (wired channels)
<<<<<<< HEAD
stalist=[STATION('COSA','AF','10','915E',('SHZ','SH1','SH2')),
        STATION('LABE','AF','10','915C',('SHZ','SHN','SHE')),
        STATION('COVA','AF','10','BAD5',('SHZ','SHN','SHE')),
        STATION('EORO','AF','10','A970',('SHZ','SHN','SHE')),
        STATION('FRAN','AF','10','A852',('SHZ','SH1','SH2')),
        STATION('GOVA','AF','10','916O',('SHZ','SHN','SHE')),
        STATION('LARB','AF','10','9C99',('SHZ','SHN','SHE')),
        STATION('MTBA','AF','10','BAEE',('SHZ','SHN','SHE')),
        STATION('MTFE','AF','10','9F9A',('SHZ','SHN','SHE')),
        STATION('MTFO','AF','10','915F',('SHZ','SHN','SHE')),
        STATION('NOLA','AF','10','9F97',('SHZ','SHN','SHE')),
        STATION('POCR','AF','10','A979',('SHZ','SH1','SH2')),
        STATION('SOLU','AF','10','BAD9',('SHZ','SHN','SHE')),
        STATION('WHAT','AF','10','A895',('SHZ','SH1','SH2')),
        STATION('REYN','AF','10','',('SHZ','SH1','SH2')),
        STATION('WHYM','AF','10','AC14',('SHZ','SHN','SHE'))]

# stalist=[STATION('ASPR','CO','10','915A',('HHZ','HH1','HH2')),
        # STATION('HUVA','CO','10','9F95',('HHZ','HH1','HH2')),
        # STATION('KING','CO','10','915D',('HHZ','HH1','HH2')),
        # STATION('MORV','CO','10','915B',('HHZ','HH1','HH2')),
        # STATION('NOBU','CO','10','AC09',('HHZ','HH1','HH2')),
        # STATION('STBA','CO','10','AC40',('HHZ','HH1','HH2')),
        # STATION('TEPE','CO','10','9F96',('HHZ','HH1','HH2'))]
=======
# stalist=[STATION('COSA','AF','10','915E',('SHZ','SH1','SH2')),
        # STATION('LABE','AF','10','915C',('SHZ','SHN','SHE')),
        # STATION('COVA','AF','10','BAD5',('SHZ','SHN','SHE')),
        # STATION('EORO','AF','10','A970',('SHZ','SHN','SHE')),
        # STATION('FRAN','AF','10','A852',('SHZ','SH1','SH2')),
        # STATION('GOVA','AF','10','916O',('SHZ','SHN','SHE')),
        # STATION('LARB','AF','10','9C99',('SHZ','SHN','SHE')),
        # STATION('MTBA','AF','10','BAEE',('SHZ','SHN','SHE')),
        # STATION('MTFE','AF','10','9F9A',('SHZ','SHN','SHE')),
        # STATION('MTFO','AF','10','915F',('SHZ','SHN','SHE')),
        # STATION('NOLA','AF','10','9F97',('SHZ','SHN','SHE')),
        # STATION('POCR','AF','10','A979',('SHZ','SH1','SH2')),
        # STATION('SOLU','AF','10','BAD9',('SHZ','SHN','SHE')),
        # STATION('WHAT','AF','10','A895',('SHZ','SH1','SH2')),
        # STATION('REYN','AF','10','',('SHZ','SH1','SH2')),
        # STATION('WHYM','AF','10','AC14',('SHZ','SHN','SHE'))]

stalist=[STATION('ASPR','CO','10','915A',('HHZ','HH1','HH2')),
        STATION('HUVA','CO','10','9F95',('HHZ','HH1','HH2')),
        STATION('KING','CO','10','915D',('HHZ','HH1','HH2')),
        STATION('MORV','CO','10','915B',('HHZ','HH1','HH2')),
        STATION('NOBU','CO','10','AC09',('HHZ','HH1','HH2')),
        STATION('STBA','CO','10','AC40',('HHZ','HH1','HH2')),
        STATION('TEPE','CO','10','9F96',('HHZ','HH1','HH2'))]
>>>>>>> 4eab03d65c72a9dc259ba578e5b556d1ab05752a

# stalist=[STATION('WMSZ','DF','10','953D',('HHZ','HHN','HHE')),
       # STATION('WPSZ','DF','10','92D1',('HHZ','HHN','HHE')),
       # STATION('WDSZ','DF','10','9337',('HHZ','HHN','HHE')),
       # STATION('WTSZ','DF','10','9793',('HHZ','HHN','HHE')),
       # STATION('WQSZ','DF','10','A895',('HHZ','HHN','HHE'))]

#######################USER EIDT SECTION, SET UP STATIONS YOU NEED TO DOWNLOAD#
geostalist=[STATION('WKZ','NZ','10','',('HHZ','HHE','HHN')),
            STATION('MLZ','NZ','10','',('HHZ','HHE','HHN')),
            STATION('JCZ','NZ','10','',('HHZ','HHE','HHN')),
            STATION('EAZ','NZ','10','',('HHZ','HHE','HHN')),
            STATION('MSZ','NZ','10','',('HHZ','HHE','HHN'))]
