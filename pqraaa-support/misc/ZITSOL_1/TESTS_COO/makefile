# this makefile is for LINUX machines only 
# 

LINKS =  -L../ -lzitsol \
	 -L/project/scicom/scicom00/SOFT/lib/Linux32 -llapack \
	 -L/project/scicom/scicom00/SOFT/lib/Linux32 -lblas 

##LINKS = -L/home/saad/scicom/lib/ -llapack \
##        -L/home/saad/scicom/lib/ -lblas

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

## ziluc  REMOVED. 

all: zarms.ex ziluk.ex zilut.ex 

zarms.ex: mainARMScoo.o   
	$(LD) $(LDFLAGS) mainARMScoo.o $(LINKS) -o zarms.ex

ziluk.ex: mainILUKcoo.o   
	$(LD) $(LDFLAGS) mainILUKcoo.o  $(LINKS) -o ziluk.ex

zilut.ex: mainILUTcoo.o   
	$(LD) $(LDFLAGS) mainILUTcoo.o $(LINKS) -o zilut.ex

##ziluc.ex: mainILUCcoo.o   
##	$(LD) $(LDFLAGS) mainILUCcoo.o $(LINKS) -o ziluc.ex

#
clean :
	rm -f *.o *.ex *~ core *.cache OUT/*
