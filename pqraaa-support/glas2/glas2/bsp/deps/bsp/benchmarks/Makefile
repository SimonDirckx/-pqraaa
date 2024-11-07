
include ../include.mk
include ../version.mk

plot_results.m: machine.info put hpput get hpget mpiput mpiget mpicollectives bsp_nprocs src/plot_results.header src/plot_results.footer
	rm -f /tmp/mcbsp_bench.log
	./put   | tee -a /tmp/mcbsp_bench.log
	./hpput | tee -a /tmp/mcbsp_bench.log
	./get   | tee -a /tmp/mcbsp_bench.log
	./hpget | tee -a /tmp/mcbsp_bench.log
	mpirun -np `./bsp_nprocs` ./mpiput | tee -a /tmp/mcbsp_bench.log
	mpirun -np `./bsp_nprocs` ./mpiget | tee -a /tmp/mcbsp_bench.log
	mpirun -np `./bsp_nprocs` ./mpicollectives | tee -a /tmp/mcbsp_bench.log
	cp src/plot_results.header ./plot_results.m
	grep ms /tmp/mcbsp_bench.log | cut -d' ' -f 4 >> ./plot_results.m
	cat src/plot_results.footer >> ./plot_results.m

../include.mk: ../Makefile
	(cd ..; make include.mk)

../lib/libmcbsp${VERSION}.a:
	(cd ..; make)

../include:
	(cd ..; make)

machine.info: generate_machine_info
	./generate_machine_info
	touch machine.info

put: src/collectives.c ../lib/libmcbsp${VERSION}.a ../include/
	${CC} ${CFLAGS} -DPRIMITIVE=0 -I../include/ -o put src/collectives.c ../lib/libmcbsp${VERSION}.a `../deplibs.sh ${CC}`

hpput: src/collectives.c ../lib/libmcbsp${VERSION}.a ../include/
	${CC} ${CFLAGS} -DPRIMITIVE=1 -I../include/ -o hpput src/collectives.c ../lib/libmcbsp${VERSION}.a `../deplibs.sh ${CC}`

get: src/collectives.c ../lib/libmcbsp${VERSION}.a ../include/
	${CC} ${CFLAGS} -DPRIMITIVE=2 -I../include/ -o get src/collectives.c ../lib/libmcbsp${VERSION}.a `../deplibs.sh ${CC}`

hpget: src/collectives.c ../lib/libmcbsp${VERSION}.a ../include/
	${CC} ${CFLAGS} -DPRIMITIVE=3 -I../include/ -o hpget src/collectives.c ../lib/libmcbsp${VERSION}.a `../deplibs.sh ${CC}`

mpiput: src/collectives.c ../lib/libmcbsp${VERSION}.a ../include/
	mpicc -std=c99 -mtune=native -march=native ${CFLAGS} -DPRIMITIVE=10 -I../include/ -o mpiput src/collectives.c

mpiget: src/collectives.c ../lib/libmcbsp${VERSION}.a ../include/
	mpicc -std=c99 -mtune=native -march=native ${CFLAGS} -DPRIMITIVE=11 -I../include/ -o mpiget src/collectives.c

mpicollectives: src/collectives.c ../lib/libmcbsp${VERSION}.a ../include/
	mpicc -std=c99 -mtune=native -march=native ${CFLAGS} -DPRIMITIVE=12 -I../include/ -o mpicollectives src/collectives.c

generate_machine_info: ../lib/libmcbsp${VERSION}.a src/generate_machine_info.cpp ../include/
	${CPP} -Wextra -O1 -g -I../include/ -o generate_machine_info src/generate_machine_info.cpp ../lib/libmcbsp${VERSION}.a `../deplibs.sh ${CPP}`

bsp_nprocs: ../lib/libmcbsp${VERSION}.a src/bsp_nprocs.c ../include/
	${CC} ${CFLAGS} -I../include/ -o bsp_nprocs src/bsp_nprocs.c ../lib/libmcbsp${VERSION}.a `../deplibs.sh ${CC}`

clean:
	rm -f machine.info generate_machine_info bsp_nprocs put hpput get hpget mpiput mpiget mpicollectives

