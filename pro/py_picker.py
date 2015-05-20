#!/usr/bin/python
"""
Python automated seismic picker implimentations for the rt2detection package.
"""

def seism_picker(picktype,args,stream):
    """
    A simple cover function for th eobspy picking modules

    :type picktype: String
    :param picktype: Eother Baer or AR for Auto-regressive
    :type args: Class
    :param args: List of appropriate arguments for the chosen pick type -
                 defined in picker_par.py
    :type Stream: obspy.Stream

    :return: p_pick, s_pick
    """
    import sys
    if picktype == 'Baer':
        from obspy.signal.trigger import pkBaer
        p_pick=[]
        phase_info=[]
        station_info=[]
        for tr in stream:
            trp_pick, trphase_info = pkBaer(tr.data, tr.stats.sampling_rate,
                                        args.tdownmax, args.tupevent, args.thr1,
                                        args.thr2, args.preset_len, args.p_dur)
                                        # 20, 60, 7.0, 12.0, 100, 100)
            trp_pick=trp_pick/tr.stats.sampling_rate
            p_pick+=[trp_pick]
            phase_info+=[trphase_info]
            station_info+=[tr.stats.station+'.'+tr.stats.channel]
        return p_pick, phase_info, station_info
    elif picktype == 'AR':
        if len(stream) != 3:
            print 'Stream must be three channels from the same station'
            sys.exit()
        from obspy.signal.trigger import arPick
        tr1=stream.select(channel='*Z')[0]
        try:
            tr2=stream.select(channel='*N')[0]
        except:
            tr2=stream.select(channel='*1')[0]
        try:
            tr3=stream.select(channel='*E')[0]
        except:
            tr3=stream.select(channel='*2')[0]
        df=tr1.stats.sampling_rate
        p_pick, s_pick = arPick(tr1.data, tr2.data, tr3.data, df,
                                args.f1, args.f2, args.lta_p, args.sta_p,
                                args.lta_s, args.sta_s, args.m_p, args.m_s,
                                args.l_p, args.l_s, args.s_pick)
                        #1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
        return p_pick, s_pick

def picker_plot(stream, picks, types, show=True, savefile=''):
    """
    Plotting fucntion for the picker routine - this will pot the waveforms and
    the picks to allow for pick veriication.

    :type stream: obspy.Stream
    :type picks: list of obspy.UTCDateTime
    :type types: list of string
    :type show: Boolean, optional
    :type savefile: String, optional
    """
    import matplotlib.pyplot as plt
    import datetime as dt
    import matplotlib.dates as mpdates
    import numpy as np
    stream.filter('highpass',freq=2.0, corners=2, zerophase=True)
    plotno=1
    fig=plt.figure(num=1, dpi=100, facecolor='w', edgecolor='k')
    fig.suptitle('Picks')
    for tr in stream:
        dates=range(0,len(tr.data))
        dates=[date/tr.stats.sampling_rate for date in dates]
        ax=fig.add_subplot(len(stream),1, plotno)
        ax.plot(dates,tr.data,'k-')
        plt.text(-0.001*len(stream[0].data),np.mean(tr.data),tr.stats.station+' '+tr.stats.channel)
        if types[plotno-1]=='P': #or types[plotno-1][1]=='P':
            ax.plot((picks[plotno-1], picks[plotno-1]), \
                    (min(tr.data)-10, max(tr.data)+10), 'r-', linewidth=2.0)
        elif types[plotno-1]=='S':
             ax.plot((picks[plotno-1], picks[plotno-1]), \
                    (min(tr.data)-10, max(tr.data)+10), 'b-', linewidth=2.0)
        plotno+=1
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
    plt.xlabel('Time [s]')
    if show:
        plt.show()
    else:
        plt.savefig(savefile)
    return

if __name__=='__main__':
    import sys, os
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print 'Requires two arguments, picker type (either AR or Baer) and the wavefile path'
        sys.exit()
    else:
        picktype=str(sys.argv[1])
        wavfile=str(sys.argv[2])
        if len(sys.argv)==4:
            station=str(sys.argv[3])
    from obspy import read
    stream=read(wavfile)
    if 'station' in locals():
        stream=stream.select(station=station)
    import sys
    sys.path.insert(0,"/home/calumch/my_programs/Built/rt2detection")
    from par import picker_par as defaults
    if picktype=='AR':
        args=defaults.AR_args
        p_pick, s_pick=seism_picker(picktype,args,stream)
        print 'P-pick at: '+str(stream[0].stats.starttime+p_pick)
        print 'S-pick at: '+str(stream[0].stats.starttime+s_pick)
        picker_plot(stream,[p_pick,s_pick,s_pick],['P','S','S'])
    elif picktype=='Baer':
        args=defaults.Baer_args
        p_pick, phase_info, station_info=seism_picker(picktype,args,stream)
        pick_types=[]
        for i in range(0,len(p_pick)):
            print 'P-pick at: '+str(stream[0].stats.starttime+p_pick[i])+\
                    ' of type: '+phase_info[i]+' on channel '+station_info[i]
            if len(phase_info[i])==4:
                pick_types+=[phase_info[i][1]]
            else:
                pick_types+=['NA']
        print (len(stream), len(p_pick), len(pick_types))
        picker_plot(stream,p_pick,pick_types)
    else:
        print 'Picktype is not AR or Baer, and it must be'
        sys.exit()
