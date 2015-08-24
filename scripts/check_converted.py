#!/usr/bin/python
"""
Script to check that all the data that are available in raw format have been
converted fully to daylong miniseed.
"""
import sys
sys.path.append('..')
from rt2detection.par import overallpars as defaults
import glob, matplotlib.pyplot as plt, datetime as dt, matplotlib.dates as mdates

rawdays=glob.glob(defaults.indir+'/*/[0-9][0-9][0-9][0-9][0-9][0-9][0-9]') #should be yyyyjjj directories
rawdays=[rawday.split('/')[-1] for rawday in rawdays]
rawdays=list(set(rawdays))
rawdays.sort()
rawdays=[dt.datetime.strptime(rawday, '%Y%j') for rawday in rawdays]
print 'There are '+str(len(rawdays))+' days of raw data'
rawstas=[] # Raw stations, will be a list of numbers per day
procstas=[] # Processed stations, will be a list of numbers per day
convstas=[] # Converted stations, will be a list of numnbers per day
f=open('Missing_das.txt','w')

for day in rawdays:
    print 'Working on day: '+day.strftime('%Y-%j')
    rawstas.append(len(glob.glob(defaults.indir+\
                                 '/*/'+day.strftime('%Y%j')+'/*/1')))
    rawdas=[sta.split('/')[-2] for sta in\
                 glob.glob(defaults.indir+'/*/'+day.strftime('%Y%j')+'/*/1')]
    conv_das=list(set([file.split('/')[-1].split('.')[-4] \
                         for file in glob.glob(defaults.outdir+'/'+\
                                               day.strftime('Y%Y/R%j.01')+\
                                               '/*')]))
    # proc_das is a list of DAS numbers that have been processed for that day
    if conv_das:
        convstas.append(len(conv_das))
        for das in rawdas:
            if not das in conv_das:
                f.write(day.strftime('%Y-%j')+' '+das+'\n')
    else:
        convstas.append(0)
        for das in rawdas:
            f.write(day.strftime('%Y-%j')+' '+das+'\n')
    procstas.append(len(list(set([file.split('/')[-1].split('.')[0] \
                                  for file in glob.glob(defaults.arcdir+'/'+\
                                                        day.strftime('Y%Y/R%j.01')+\
                                                        '/*.AF.*')]))))

f.close()
ax=plt.subplot(111)
plt.plot(rawdays, rawstas, 'k', linewidth=2, label='Raw data')
plt.plot(rawdays, convstas, 'b', linewidth=2, label='Converted data')
plt.plot(rawdays, procstas, 'r', linewidth=2, label='Processed data')
dateFormatter = mdates.DateFormatter('%Y-%m-%d')
ax.xaxis.set_major_formatter(dateFormatter)
plt.xlabel('Date')
plt.ylabel('Number of channels')
plt.legend()
plt.show()
