"""
rt2detection is a simple set of functions to convert reftek data to multiplexed
miniseed data for seisan integration.  This routine will also run a simple
sta/lta routine through the data to detect events.  It can then run the
filterpicker routine of Anthony Lomax over detected events and locate these
events using hypocentre.  In the end you may be able to generate a reasonable
catalogue just through these routines, but care must be taken when selecting
parameters and remember that automatic picks need checking - filterpicker
is a fast but innacurate picker.  Be warned.

"""

import sys
sys.path.insert(0,"home/calumch/my_programs/Building/rt2detection")

__all__ = ['par', 'pro']

if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)
