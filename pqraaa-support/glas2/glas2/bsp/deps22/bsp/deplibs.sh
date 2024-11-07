#!/bin/bash
#convenience script, used by Makefile to enable easy cross-platform linking (Linux <> OS X):
#returns the library dependencies of MulticoreBSP for C on this platform
#uses a trick from user mhawke at stackoverflow.com
if [ $# -eq 0 ]
then
	CC=gcc
else
	CC=$1
fi

echo "int main(){}" | ${CC} -o /dev/null -x c - -lrt 2>/dev/null && echo "-lpthread -lrt" || echo "-lpthread"
