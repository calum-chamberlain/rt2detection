#!/usr/bin/python
"""
Codes to extract amplitude information from s-file and cnvert these to a
previously dervied local magnitude scale.  Originally written to generate
magnitudes for the DFDP Alpine Fault drilling project in New Zealand by
Calum Chamberlain (VUW) using the magnitude constasts derived by Boese et al.
(2012).
"""

# import parameter file
import numpy as np

def dist_calc(loc1, loc2):
    """
    Function to calcualte the distance in km between two points, uses the flat
    Earth approximation

    :type loc1: Tuple
    :param loc1: Tuple of lat, lon, depth (in decimal degrees and km)
    :type loc2: Tuple
    """
    R=6371.009  # Radius of the Earth in km
    dlat=np.radians(abs(loc1[0]-loc2[0]))
    dlong=np.radians(abs(loc1[1]-loc2[1]))
    ddepth=abs(loc1[2]-loc2[2])
    mean_lat=np.radians((loc1[0]+loc2[0])/2)
    dist=R*np.sqrt(dlat**2+(np.cos(mean_lat)*dlong)**2)
    dist=np.sqrt(dist**2+ddepth**2)
    return dist

def mag_conv(deltaamp, distance, sta_cor, f_dependent=False, period=999.9):
    """
    Function to convert from peak-to-trough amplitude to magnitude for a single
    seismic station at a certain hypocentral distance (distance) and with a set
    station correction term.  Will look for parameters for conversion in a file
    mag_conv_par,py in the ../par directory (relative to this code).

    Certain assumptions are implicit in this calculation:
        *   peak-to-peak amplitude represents the maximum peak-to-trough
            amplitude of a body wave phase on a Wood-Anderson converted
            seismogram and is measured in nm.  This will be halved to give an
            approximate peak amplitude as specified in the definition of the
            Wood Anderson local magnitude scale.

    :type deltaamp: float
    :param deltaamp: peak-to-peak amplitude in nm
    :type distance: float
    :param distance: hypocentral distance in km
    :type sta_corr: float
    :param sta_corr: Station correction term in magnitude units.
    :type f_dependent: bool
    :param f_dependent: option, defaults to False
    :type period: float
    :param period: optional, required if f_dependent=True

    :return: Magnitude, float
    """
    from par import mag_conv_par as mag_par

    # Perform a basic check to confirm that the local magnitude scale applies
    if distance > mag_par.maxdist:
        raise ValueError('Earthquake beyond magnitude scale')

    if f_dependent:
        if period==999.9:
            raise ValueError('You have set frequency dependence to true, but not given a period')
        gamma=mag_par.gamma*period
    else:
        gamma=mag_par.gamma

    Ml=(np.log10(0.5*deltaamp))+\
            (mag_par.alpha*np.log10(distance))+\
            (0.4343*gamma*distance)-\
            sta_cor
    return Ml

def event_magnitude(sfile):
    """
    Function to generate the magnitude of a single event from a seisan s-file
    with amplitude picks in it

    :type sfile: String
    :param sfile: Nordic type s-file name, with full path

    :return: Local magnitude, standard deviation of magnitude
    """
    from pro import Sfile_util
    from par import mag_conv_par as mag_par
    # Check that the s-file exists
    import glob
    if not glob.glob(sfile):
        raise NameError('Sfile does not exist: '+sfile)

    picks=Sfile_util.readpicks(sfile)
    Mag_out=[]
    for pick in picks:
        if pick.phase=='IAML':
            if 'sta_cor' in locals():
                del sta_cor
            for station in mag_par.station_corrections:
                if pick.station==station[0]:
                    sta_cor=station[1]
            if not 'sta_cor' in locals():
                sta_cor=1.0
                # print '\nStation correction not found for station '+pick.station
            try:
                Magnitude=mag_conv(pick.amplitude, pick.distance, sta_cor,\
                                  mag_par.frequency_dependent,\
                                  pick.peri)
                if not np.isnan(Magnitude):
                    Mag_out.append(Magnitude)
            except (ValueError):
                # print 'Either earthquake too far away, or frequency is not set'
                pass

    Mag_std=np.std(Mag_out)
    Mag_out=np.mean(Mag_out) # Take the mean magnitude
    return Mag_out, Mag_std

def recalc_database(path, plot=True):
    """
    Overarching code to recalculate and plot all the magnitudes for a given
    database in seisan.  Must be a databse of S-files in a SEISAN REA structure

    :type path: String
    :param path: Path to the top of the rea tree (above the year directories)

    :return: Event info, list of tuples (Mag_in, Mag_out, Date, Location)
    """
    import glob, sys
    from par import mag_conv_par as mag_par
    from pro import Sfile_util
    if not glob.glob(path):
        raise NameError('Path does not exist '+path)

    sfilelist=glob.glob(path+'/*/*/*.S*')
    Mag_in=[]
    Mag_out=[]
    Date=[]
    Event_info=[]
    for sfile in sfilelist:
        sys.stdout.write('Working on sfile: '+sfile+'\r')
        sys.stdout.flush()
        Date.append(Sfile_util.readheader(sfile).time)
        if not np.isnan(Sfile_util.readheader(sfile).Mag_1):
            Mag_in.append(Sfile_util.readheader(sfile).Mag_1)
        Magnitude=(event_magnitude(sfile)[0])
        if not np.isnan(Magnitude):
            Mag_out.append(Magnitude)
        Event_info.append([Magnitude, Sfile_util.readheader(sfile).Mag_1, \
                Sfile_util.readheader(sfile).time, \
                (Sfile_util.readheader(sfile).latitude, \
                Sfile_util.readheader(sfile).longitude, \
                Sfile_util.readheader(sfile).depth)])
    if plot:
        import matplotlib.pyplot as plt
        try:
            # Plot histogram
            n, bins, patches=plt.hist(Mag_in,len(Mag_in)/5, facecolor='Black', \
                                      alpha=0.5, label='Previous, n='+str(len(Mag_in)))
            n, bins, patches=plt.hist(Mag_out,len(Mag_out)/5, facecolor='Red', \
                                      alpha=0.75, label='Recalculated, n='+str(len(Mag_out)))
            plt.legend()
            plt.ylabel('Number of events')
            plt.xlabel('Local Magnitude $M_L$')
            plt.show()
            Mag_out=np.sort(Mag_out)
            cdf=np.arange(len(Mag_out))/float(len(Mag_out)) # normalized, useful in a mo
            cdf=((cdf*-1.0)+1.0)*len(Mag_out)
            plt.plot(Mag_out,np.log10(cdf), 'r', linewidth=2.0, label='Recalculated')
            Mag_in=np.sort(Mag_in)
            cdf=np.arange(len(Mag_in))/float(len(Mag_in)) # normalized, useful in a mo
            cdf=((cdf*-1.0)+1.0)*len(Mag_in)
            plt.plot(Mag_in,np.log10(cdf), 'k', linewidth=2.0, label='Previous')
            plt.legend()
            plt.ylabel('$Log_{10}$ of cumulative density')
            plt.xlabel('Local Magnitude$M_L$')
            plt.show()
            return Event_info
        except (AttributeError):
            print 'Error plotting'
            return Event_info
    else:
        return Event_info


def plot_dist_mag(path, origin):
    """
    Function to recalculate the magnitudes for a given database and plot the
    magnitude against distance.

    :type path: Str
    :param path: Database to convert
    :type origin: Tuple
    :param origin: Lat, Long and Depth of origin
    """
    Event_info=recalc_database(path, False)
    mag_in=[]
    mag_new=[]
    dist_in=[]
    dist_out=[]
    for event in Event_info:
        if not np.isnan(event[0]):
            mag_new.append(event[0])
            dist_out.append(dist_calc(origin, event[3]))
        if not np.isnan(event[1]):
            mag_in.append(event[1])
            dist_in.append(dist_calc(origin, event[3]))
    import matplotlib.pyplot as plt
    plt.scatter(dist_out,mag_new, c='red', label='Recaclulated')
    plt.scatter(dist_in,mag_in,c='black', alpha=0.3, label='Original')
    plt.xlabel('Distance in km')
    plt.ylabel('Magnitude')
    plt.xlim([0,350])
    plt.title('Magnitude with distance from '+str(origin[0])+', '+str(origin[1]))
    # plt.legend()
    plt.show()
    return

def plot_dist_RMS(path, origin):
    """
    Function to recalculate the magnitudes for a given database and plot the
    magnitude against distance.

    :type path: Str
    :param path: Database to convert
    :type origin: Tuple
    :param origin: Lat, Long and Depth of origin
    """
    import glob
    from pro import Sfile_util
    sfilelist=glob.glob(path+'/*/*/*.S*')
    RMS=[]
    dist=[]
    for sfile in sfilelist:
        header=Sfile_util.readheader(sfile)
        RMS.append(header.t_RMS)
        dist.append(dist_calc(origin, (header.latitude, header.longitude, header.depth)))
    import matplotlib.pyplot as plt
    # plt.semilogy(dist,RMS, c='red', marker='o', ls='None')
    plt.plot(dist,RMS, c='red', marker='o', ls='None')
    plt.xlabel('Distance in km')
    plt.ylabel('RMS (s)')
    plt.xlim([0,350])
    plt.ylim([0,10])
    plt.title('RMS with distance from '+str(origin[0])+', '+str(origin[1]))
    # plt.legend()
    plt.show()
    return

def plot_residuals(path):
    """
    Function to read in S-files and make a histogram of the pick residuals for
    each station

    :type path: Str
    """
    import glob, sys
    import matplotlib.pyplot as plt
    from pro import Sfile_util
    import numpy as np
    sfilelist=glob.glob(path+'/*/*/*.S*')
    if 'stachan_list' in locals():
        del stachan_list
    for sfile in sfilelist:
        picks=Sfile_util.readpicks(sfile)
        for pick in picks:
            match=0
            if not 'stachan_list' in locals():
                if abs(pick.timeres) <=4.0 and pick.channel != 'HT':
                    if pick.phase=='P':
                        stachan_list=[[pick.station, pick.channel, [pick.timeres],[]]]
                    elif pick.phase=='S':
                        stachan_list=[[pick.station, pick.channel, [], [pick.timeres]]]
            else:
                for stachan in stachan_list:
                    if pick.station == stachan[0] and \
                            pick.channel == stachan[1] and not np.isnan(pick.timeres)\
                            and abs(pick.timeres) <= 4.0:
                        if pick.phase=='P':
                            stachan[2].append(pick.timeres)
                        elif pick.phase=='S':
                            stachan[3].append(pick.timeres)
                        match=1
                if match == 0 and not np.isnan(pick.timeres)\
                        and abs(pick.timeres) <= 4.0 and pick.channel != 'HT':
                    print 'Found picks for: '+pick.station+' '+pick.channel
                    if pick.phase=='P':
                        stachan_list.append([pick.station, pick.channel, [pick.timeres], []])
                    elif pick.phase=='S':
                        stachan_list.append([pick.station, pick.channel, [], [pick.timeres]])
    # Print out some useful stats
    ppicks=0
    spicks=0
    presidual=0
    sresidual=0
    for stachan in stachan_list:
        ppicks+=len(stachan[2])
        spicks+=len(stachan[3])
        presidual+=sum(stachan[2])
        sresidual+=sum(stachan[3])
    print 'Total P-picks: '+str(ppicks)
    print 'Total S-picks: '+str(spicks)
    print 'P RMS mean: '+str(presidual/ppicks)
    print 'S RMS mean: '+str(sresidual/spicks)
    # Plot the results
    i=0
    # Get unique list of stations, make one plot for each station
    stations=[]
    for stachan in stachan_list:
        stations+=[stachan[0]]
    stations=list(set(stations))
    # Concatenate all the picks for each station
    stachan_list.sort()
    for stachan in stachan_list:
        if not 'sta_list' in locals():
            sta_list=[[stachan[0], 'all', stachan[2], stachan[3]]]
            station=stachan[0]
            i=0
        else:
            if station==stachan[0]:
                sta_list[i][3]+=stachan[3]
                sta_list[i][2]+=stachan[2]
            else:
                sta_list.append(stachan)
                i+=1
                station=stachan[0]
    fig, axes = plt.subplots((len(sta_list)), 1, sharex=True)#, sharey=True)
    print 'I have picks for '+str(len(sta_list))+' stations'
    axes=axes.ravel()
    i=0
    for stachan in sta_list:
        print 'Plotting for station: '+stachan[0]
        if len(stachan[2]) != 0:
            n, bins, patches=axes[i].hist(stachan[2], bins=np.arange(-4.0, 4.0, 0.025)\
                                          , facecolor='Black', alpha=0.5)
            axes[i].text(0.85, 0.8, r'$\ P:\ \mu='+str(np.mean(stachan[2]))[0:4]+\
                         ',\ \sigma='+str(np.std(stachan[2]))[0:4]+\
                         ',\ n='+str(len(stachan[2]))+'$',\
                         horizontalalignment='center', verticalalignment='center',
                         transform=axes[i].transAxes)
        if len(stachan[3]) != 0:
            n, bins, patches=axes[i].hist(stachan[3],bins=np.arange(-4.0, 4.0, 0.025)\
                                          , facecolor='Red', alpha=0.75)
            axes[i].text(0.15, 0.8, r'$\ S:\ \mu='+str(np.mean(stachan[3]))[0:4]+\
                         ',\ \sigma='+str(np.std(stachan[3]))[0:4]+\
                         ',\ n='+str(len(stachan[3]))+'$',\
                         horizontalalignment='center', verticalalignment='center',
                         transform=axes[i].transAxes, color='Red')
        axes[i].set_ylabel(stachan[0])
        # axes[i].yaxis.set_label_position("right")
        axes[i].yaxis.tick_right()
        axes[i].locator_params(axis='y', nbins=2)
        i+=1
    axes[i-1].set_xlabel('RMS residual (s)')
    plt.xlim(-2.0, 2.0)
    # plt.ylim(0,40)
    fig.subplots_adjust(hspace=0.25)
    fig.subplots_adjust(wspace=0)
    fig.text(0.94, 0.5, 'Number of picks', ha='center', va='center', rotation=270)
    plt.show()
    # plt.savefig('residuals.eps')
    return stachan_list

def summary_table(path, output='csv'):
    """
    Function to generate a summary table of earthquake information. Can output
    either as a .csv or a .tex file.

    :type path: str
    :param path: Database directory
    :type output: str
    :param output: Either 'csv' or 'tex'
    """
    # Check that the output is correct
    if not output in ['csv', 'tex']:
        print output
        raise ValueError('Output format not recognised')
    import glob
    from pro import Sfile_util
    sfilelist=glob.glob(path+'/*/*/*.S*')
    f=open('Summary.'+output, 'w')
    # Write header
    if output == 'csv':
        f.write('Date, Origin time (UTC), Latitude (deg), Longitude (deg), Depth (km), Magnitude(seisan), Magnitude(local)\n')
    elif output == 'tex':
        f.write('\begin{table}{c c c c c c c}\n')
        f.write('\textbf{Date} & \textbf{Origin time (UTC)} &  \textbf{Latitude (deg)}'+\
                '\textbf{Longitude (deg)} & \textbf{Depth (km)} & \textbf{Magnitude(seisan)}'+\
                '\textbf{Magnitude (local)}\\ \n')
        f.write('\hline\n')
    # Write contents
    for sfile in sfilelist:
        EQ_info=Sfile_util.readheader(sfile)
        Mag_out, Mag_std = event_magnitude(sfile)
        if output == 'csv':
            f.write(str(EQ_info.time.year)+'/'+str(EQ_info.time.month).zfill(2)+'/'+\
                    str(EQ_info.time.day).zfill(2)+', '+\
                    str(EQ_info.time.hour).zfill(2)+':'+\
                    str(EQ_info.time.minute).zfill(2)+':'+\
                    str(EQ_info.time.second).zfill(2)+'.'+\
                    str(EQ_info.time.microsecond).zfill(2)+', '+\
                    str(EQ_info.latitude)+', '+str(EQ_info.longitude)+', '+\
                    str(EQ_info.depth)+','+str(EQ_info.Mag_1)+','+\
                    str(Mag_out)+'\n')
        elif output == 'tex':
            f.write(str(EQ_info.time.year)+'/'+str(EQ_info.time.month).zfill(2)+'/'+\
                    str(EQ_info.time.day).zfill(2)+' & '+\
                    str(EQ_info.time.hour).zfill(2)+':'+\
                    str(EQ_info.time.minute).zfill(2)+':'+\
                    str(EQ_info.time.second).zfill(2)+'.'+\
                    str(EQ_info.time.microsecond).zfill(2)+' & '+\
                    str(EQ_info.latitude)+' & '+str(EQ_info.longitude)+' & '+\
                    str(EQ_info.depth)+' & '+str(EQ_info.Mag_1)+' & '+\
                    str(Mag_out)+'\\ \n')
    # Write end of table for latex
    if output == 'tex':
        f.write('\end{table}')
    f.close()
    print 'Written summary file: Summary.'+output
    return

if __name__ == '__main__':
    import sys
    if len(sys.argv)!=2:
        print 'Requires one argument, the path to the s-file database'
        sys.exit()
    else:
        sys.path.insert(0,"/home/calumch/my_programs/Building/rt2detection")
        from par import mag_conv_par as mag_par
        from pro import Sfile_util
        recalc_database(str(sys.argv[1]))
        print '\n'
