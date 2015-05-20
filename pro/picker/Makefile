
DIRS = libmseed src

all rtquake clean static install gcc ::
	@for d in $(DIRS) ; do \
	    echo "Running $(MAKE) $@ in $$d" ; \
	    if [ -f $$d/Makefile -o -f $$d/makefile ] ; \
	        then ( cd $$d && $(MAKE) $@ ) ; \
	    elif [ -d $$d ] ; \
	        then ( echo "ERROR: no Makefile/makefile in $$d for $(CC)" ) ; \
	    fi ; \
	done

