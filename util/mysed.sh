#!/bin/bash
for file in `ls *.c *.h` ;
do
	echo $file
	sed s/"$1"/"$2"/ $file > abc
	cp abc $file
done
