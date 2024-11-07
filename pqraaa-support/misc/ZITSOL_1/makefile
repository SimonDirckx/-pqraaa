# this makefile is for LINUX machines only 
#  ## iluc removed 
OBJ	=  zfgmres.o ziluk.o zilut.o zarms2.o auxill.o

LIB  = ./LIB/zPQ.o ./LIB/zpiluNEW.o ./LIB/zindsetC.o \
       ./LIB/zsets.o ./LIB/ztools.o ./LIB/systimer.o ./LIB/zmisc.o  \
       ./LIB/zMatOps.o ./LIB/zilutpC.o ./LIB/givens.o

AR = ar -rcv

#
FC      =  f77
FCFLAGS =  -c -g -Wall
CC      =  cc
CCFLAGS =  -c -DLINUX -Wall -O3
LD      =  f77 
LDFLAGS = 
#
# clear list of default suffixes, and declare default suffixes
.SUFFIXES:
.SUFFIXES: .f .c .o
# default rule to make .o files from .f files
.f.o  : ;       $(FC) $(FCFLAGS) $*.f -o $*.o
.c.o  : ;       $(CC) $(CCFLAGS) $*.c -o $*.o
#

lib libzitsol.a: $(OBJ) $(LIB)
	$(AR) libzitsol.a  $(OBJ) $(LIB)
##       ranlib libzitsol.a
	     
#
clean :
	rm -f ${OBJ} *.o *.ex *~ core ${LIB} LIB/*~ OUT/*

cleanall :
	rm -f ${OBJ} *.o *.ex *.a *.cache *~ core ${LIB} \
        LIB/*~ LIB/*.cache OUT/* TESTS_COO/*.o TESTS_COO/*.ex \
        TESTS_HB/*.o TESTS_HB/*.ex 

