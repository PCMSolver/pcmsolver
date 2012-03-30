#!/bin/bash
list="WEM
WEMPCG
WEMPGMRES
WEMRHS
basis
compression
constants
cubature
data
dwt
energy
gauss_legendre
gauss_square
integrate
interpolate
intlin1
intlin2
intlin3
intlin4
intvector
kern
mask
phi
postproc
precond
quadrature
read_points
sparse
sparse2
topology
trafos
vector2
vector3
volume"

echo $list

for file in $list ;
do
	UPNAME=`echo $file | tr '[:lower:]' '[:upper:]'`
	first="#ifndef $UPNAME" 
	second="#define $UPNAME"
	last="#endif"
	echo $first > tempfile.txt
	echo $second >> tempfile.txt
	cat $file.h >> tempfile.txt
	echo $last >> tempfile.txt
	cp tempfile.txt $file.h
done

