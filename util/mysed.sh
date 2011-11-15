#!/bin/bash
for file in `ls *.c *.h` ;
do
	echo $file
	sed s/"malloc.h"/"stdlib.h"/ $file > abc
	cp abc $file
done
