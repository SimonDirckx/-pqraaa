
include ../../include.mk

INCLUDE_FLAGS = -I../include/ -I../../include/
LFLAGS = ../lib/libmultibsp${VERSION}.a ../../lib/libmcbsp${VERSION}.a

CPPFLAGS := ${CPPFLAGS} ${INCLUDE_FLAGS}

all: ../include ../../include multibsp-ip multibsp-fft

multibsp-ip: multibsp-ip.opp11
	${CPP11} ${CPPFLAGS} -o $@ $^ ${LFLAGS} `../../deplibs.sh ${CPP11}`

multibsp-fft: multibsp-fft.opp11
	${CPP11} ${CPPFLAGS} -o $@ $^ ${LFLAGS} `../../deplibs.sh ${CPP11}` -lfftw3

clean:
	rm -f multibsp-ip.opp11 multibsp-ip multibsp-fft.opp11 multibsp-fft

../include:
	(cd ..; make)

../../include:
	(cd ../..; make)

