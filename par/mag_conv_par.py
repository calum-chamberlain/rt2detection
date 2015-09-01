#!/usr/bin/python
"""
Simple outline of parameters for the conversion of amplitudes into magnitudes

Parameters herein are those found by Boese et al. (2012) for the Southern Alps
of New Zealand.
"""
station_corrections=[('MTFO',1.65),
                     ('WHYM',1.6),
                     ('LABE',1.45),
                     ('COVA',1.4),
                     ('COSA',1.4),
                     ('WHAT',0.86),
                     ('GOVA',1.25),
                     ('FRAN',1),
                     ('EORO',1.6),
                     ('POCR',1.12),
                     ('FOZ',1.5744),
                     ('DCZ',2.5517),
                     ('GCSZ',1.5434),
                     ('JCZ',1.7319),
                     ('LBZ',1.5777),
                     ('RPZ',1.4243),
                     ('WDSZ',1.7786),
                     ('WHFS',2.0188),
                     ('WMSZ',1.4374),
                     ('WPSZ',1.556),
                     ('WQSZ',1.3807),
                     ('WTSZ',1.709),
                     ('WVZ',1.6275)]
station_corrections=[(sta[0],sta[1]-0.5) for sta in station_corrections]
# Need other station correction terms for WQSZ, WTSZ, WDSZ etc.
alpha=1.0
frequency_dependent=True
# gamma=0.00189               # Anelastic attenuation parameter in either \km or
gamma=0.00237               # Anelastic attenuation parameter in either \km or
                            # frequency_dependent=False, or in s/km if True.
maxdist=100.0                # Maximum distance for the magnitude scale in ke
