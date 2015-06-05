
# Build environment can be configured the following
# environment variables:
#   CC : Specify the C compiler to use
#   CFLAGS : Specify compiler options to use

# Options specific for GCC
#GCC = gcc
#GCCFLAGS = -O2 -Wall -I../libmseed
BINDIR=.
CC = gcc
# CC = /usr/local/bin/gcc

# Options specific for GCC
OBJS_PICK = rdtrigL.c FilterPicker5.c FilterPicker5_Memory.c PickData.c

all :
	$(CC) -O -I.. -o rdtrigL ${OBJS_PICK} -lm -L../libmseed -lmseed -I../libmseed

rtquake:
	$(CC) -D RTQUAKE -O -I.. -o rdtrigL ${OBJS_PICK} -lgd -lm -L../libmseed -lmseed -I../libmseed

install :
	mv rdtrigL ../../bin/rtpick

clean :
	rm -f *.o *.a ${BINDIR}/rdtrigL
