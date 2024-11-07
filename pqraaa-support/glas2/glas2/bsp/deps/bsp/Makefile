include include.mk
include version.mk

ifeq ($(wildcard machine.info),)
	COPY_MACHINFO = true
else
	COPY_MACHINFO = false
endif

library: ${OBJECTS} mcbsp.opp
	mkdir include || true
	mkdir lib || true
	cp mcbsp.h include/
	cp mcbsp.hpp include/
	cp mcbsp-affinity.h include/
	cp mcbsp-templates.hpp include/
	${AR} rcs lib/libmcbsp${VERSION}.a ${OBJECTS} mcbsp.opp

compat:
	rm -f ${OBJECTS}
	make CFLAGS="${CFLAGS} -DMCBSP_COMPATIBILITY_MODE" ${OBJECTS} tests
	mkdir include || true
	mkdir lib || true
	cp mcbsp.h include/
	cp mcbsp-affinity.h include/
	${AR} rcs lib/compat-libmcbsp${VERSION}.a ${OBJECTS}
	rm -f ${OBJECTS}
	rm -f ${TESTOBJ} tests/test_abort.out

debug:
	rm -f ${OBJECTS}
	make OPT="${DEBUGFLAGS} -DMCBSP_USE_MUTEXES" ${OBJECTS} tests
	mkdir include || true
	mkdir lib || true
	cp mcbsp.h include/
	cp mcbsp-affinity.h include/
	${AR} rcs lib/libmcbsp${VERSION}-debug.a ${OBJECTS}
	rm -f ${OBJECTS} ${TESTOBJ} tests/test_abort.out
	make OPT="${DEBUGFLAGS}" CFLAGS="${CFLAGS} -DMCBSP_USE_MUTEXES -DMCBSP_COMPATIBILITY_MODE" ${OBJECTS} tests
	${AR} rcs lib/compat-libmcbsp${VERSION}-debug.a ${OBJECTS}
	rm -f ${OBJECTS} ${TESTOBJ} tests/test_abort.out

include.mk: include.default
	cp include.default include.mk

clean:
	rm -f ${OBJECTS} ${TESTOBJ} mcbsp.opp tests/test_abort.out
	rm -f test_internal    test_util    test_bsp    test_abort    tests/test_abort.out    test_spmd
	rm -rf include
	rm -rf lib

tests: mcinternal.o mcutil.o mcbsp.o ${TESTOBJ}
	${CC} ${CFLAGS} -o test_internal tests/internal.o mcinternal.o mcutil.o mcbsp.o ${LFLAGS}
	${CC} ${CFLAGS} -o test_util     tests/util.o     mcinternal.o mcutil.o mcbsp.o ${LFLAGS}
	${CC} ${CFLAGS} -o test_bsp      tests/bsp.o      mcinternal.o mcutil.o mcbsp.o ${LFLAGS}
	${CC} ${CFLAGS} -o test_abort	 tests/abort.o	  mcinternal.o mcutil.o mcbsp.o ${LFLAGS}
	${CC} ${CFLAGS} -o test_spmd     tests/spmd.o     mcinternal.o mcutil.o mcbsp.o ${LFLAGS}
	rm -f tests/internal.o tests/util.o tests/bsp.o tests/abort.o tests/test_abort.out tests/spmd.o
	if ${COPY_MACHINFO}; then cp tests/machine.info ./; fi
	./test_internal
	./test_util > /dev/null; if [ "$$?" -ne "0" ]; then echo "SUCCESS"; else echo "FAILURE"; fi
	./test_bsp
	./test_abort 2> tests/test_abort.out; diff tests/test_abort.out tests/test_abort.chk; if [ $$? -eq 0 ]; then echo "SUCCESS"; else echo "FAILURE"; fi
	./test_spmd
	if ${COPY_MACHINFO}; then rm -f machine.info; fi

