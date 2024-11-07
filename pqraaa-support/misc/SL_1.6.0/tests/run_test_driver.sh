#!/bin/bash
./test_driver &> test.out
SCHEMES=`grep -A 1000 Possible test.out | grep -B 1000 in-place | grep -v Possible | grep : | cut -d':' -f1`
if [ ! -f hardware.info ]; then
	read -p "How many cores can I run on? (Please input a positive integer) " CORES
	if [[ ${CORES} =~ ^[0-9]+$ ]]; then
		echo "${CORES}" > hardware.info
	else
		echo "Invalid input: ${CORES}. Will set number of available cores to 1...";
		echo "1" > hardware.info
	fi
fi
for i in ${SCHEMES}; do
	echo Scheme `./test_driver | grep "^[ ]*$i:"`
	./test_driver tests/s3dkt3m2/s3dkt3m2.mtx $i 1 3 1 | grep MSE
	echo
done
rm -f test.out

