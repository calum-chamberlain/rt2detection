#!/usr/bin/python

#-------------------------------------------------------------------
#   Purpose: Convenience imports for rt2detection
#   Author:  Calum J. Chamberlain
#-------------------------------------------------------------------

"""
__init__.py file for the pro subpackage containing the main routines
recquired for the package

rt2detection.pro - programs called by rt2detection
"""



__all__ = ['mag_conv']

if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)
