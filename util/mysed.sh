#!/bin/bash
for file in `ls *.c *.h` ;
do
	echo $file
	sed s/"$1"/"$1_pwl"/ $file > abc
	cp abc $file
done
