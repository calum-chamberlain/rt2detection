#!/usr/bin/python
"""
Python imput arguments for the python auto-pickers
"""
# Arguments for the Baer picker
class BAER_ARGS:
    """Input arguments for the Baer picker"""
    def __init__(self, tdownmax, tupevent, thr1, thr2,
                 preset_len, p_dur):
        self.tdownmax=tdownmax
        self.tupevent=tupevent
        self.thr1=thr1
        self.thr2=thr2
        self.preset_len=preset_len
        self.p_dur=p_dur

class AR_ARGS:
    """Input arguments for the AR picker"""
    def __init__(self, f1, f2, lta_p, sta_p, lta_s, sta_s,
                m_p, m_s, l_p, l_s, s_pick):
        self.f1=f1 # Lower bandpass frequency
        self.f2=f2 # upper bandpass frequency
        self.lta_p=lta_p #length of LTA for p arrival in seconds
        self.sta_p=sta_p # length of STA for p arrival in seconds
        self.lta_s=lta_s # length of LTA for s arrival in seconds
        self.sta_s=sta_s # length of STA for s arrival in seconds
        self.m_p=m_p # number of AR coefficients for p arrival
        self.m_s=m_s # number of AR coefficients for s arrival
        self.l_p=l_p # length of variance window for p arrival in seconds
        self.l_s=l_s # length of variance window for s arrival in seconds
        self.s_pick=s_pick # if True also pick S, otherwise only P

#################### DEFINE DEAFULT VALUES####################################
Baer_args=BAER_ARGS(20, 60, 7.0, 12.0, 100, 100)

AR_args=AR_ARGS(5.0, 30.0, 1.0, 0.1, 2.5, 0.2, 2, 8, 0.1, 0.2, False)
