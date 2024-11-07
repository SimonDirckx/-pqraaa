
include ../include.mk
include ../version.mk

all: library examples tests
	./test_sum

library: CPPFLAGS += -I../
library: multibsp.opp11
	mkdir include || true
	mkdir lib || true
	cp multibsp.hpp include/
	${AR} rcs lib/libmultibsp${VERSION}.a multibsp.opp11

examples: library multibsp_id multibsp_ip multibsp_ip_templates multibsp_fft

tests: test_sum

multibsp_id: CPPFLAGS += -I./include/ -I../include/
multibsp_id: examples/multibsp-id.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a
	${CPP11} -o multibsp_id examples/multibsp-id.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a `../deplibs.sh ${CPP11}`

multibsp_ip: CPPFLAGS += -I./include/ -I../include/
multibsp_ip: examples/multibsp-ip.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a
	${CPP11} -o multibsp_ip examples/multibsp-ip.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a `../deplibs.sh ${CPP11}`

multibsp_ip_templates: CPPFLAGS += -I./include/ -I../include/ -I./templates/
multibsp_ip_templates: examples/multibsp-ip-templates.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a
	${CPP11} -o multibsp_ip_templates examples/multibsp-ip-templates.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a `../deplibs.sh ${CPP11}`

multibsp_fft: CPPFLAGS += -I./include/ -I../include/ -I./templates/
multibsp_fft: examples/multibsp-fft.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a
	${CPP11} -o multibsp_fft examples/multibsp-fft.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a `../deplibs.sh ${CPP11}` -lfftw3 -lnuma

test_sum: CPPFLAGS += -I./include/ -I../include/ -I./templates/
test_sum: test/sum.opp11 test/hardcoded-addition.opp11 test/templated-addition.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a
	${CPP11} -o test_sum test/sum.opp11 test/hardcoded-addition.opp11 test/templated-addition.opp11 ../lib/libmcbsp${VERSION}.a ./lib/libmultibsp${VERSION}.a `../deplibs.sh ${CPP11}`

assembly: CPPFLAGS += -I./templates/ -I../include/
assembly: test/hardcoded-addition.cpp test/templated-addition.cpp
	${CPP11} ${CPPFLAGS} -S -o test/hardcoded-addition.s test/hardcoded-addition.cpp
	${CPP11} ${CPPFLAGS} -S -o test/templated-addition.s test/templated-addition.cpp

clean:
	rm -rf include
	rm -rf lib
	rm -f multibsp.opp11 examples/multibsp-ip.opp11 examples/multibsp-ip-templates.opp11 test/sum.opp11 test_sum examples/multibsp-fft.opp11 examples/multibsp-id.opp11

