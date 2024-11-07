
include ../include.mk

CPPFLAGS:=${CPPFLAGS} -I../

library: multibsp.opp11
	mkdir include || true
	mkdir lib || true
	cp initialiser.hpp include/
	cp destructor.hpp include/
	cp initialisable.hpp include/
	cp array.hpp include/
	cp distributed-vector.hpp include/
	cp replicated-vector.hpp include/
	cp leaf-data.hpp include/
	cp multibsp.hpp include/
	${AR} rcs lib/libmultibsp${VERSION}.a multibsp.opp11

clean:
	rm -rf include
	rm -rf lib
	rm -f multibsp.opp11

