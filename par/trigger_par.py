#!/usr/bin/python

"""
Declaration of variables to be used as defaults in the STALTA.py
routine.
"""
# Setting for detection routines
stalen=30   # STA length in seconds
ltalen=700  # LTA length in seconds
trigon=3.0  # sta/lta ratio to turn trigger on
trigoff=2.0 # sta/lta ratio to turn trigger off
crat=0.8    # Ratio parameter for carl sta/lta routine
cquiet=0.8  # Quiet parameter for carl sta/lta routine
netsum=3    # Stations to detect over
netwin=30   # Window in seconds for a coincidence trigger to be flagged
comp='Z'    # Component to detect on, Z,N,E or all
highcut=20  # High cut for banpdass filter in Hz
lowcut=5    # Low cut for bandpass filter in Hz
precut=90   # Time to cut before the trigger in seconds
postcut=150 # Time to cut after the trigger in seconds
<<<<<<< HEAD
net='SAMTR' # Network code, must be 5 characters, if less fill with underscores
=======
net='COSAW' # Network code, must be 5 characters, if less fill with underscores
>>>>>>> 4eab03d65c72a9dc259ba578e5b556d1ab05752a
trigout='Y' # Set to 'Y' to output triggered files, will be output locally
