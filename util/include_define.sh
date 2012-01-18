#!/bin/bash
list="BEMInterface 
GreensFunctionInterface 
WEM
WEMPCG
WEMPGMRES
WEMRHS
basis
compression
constants
cubature
dalton_wem
data
dwt
energy
gamma
gauss_legendre
gauss_square
init_points
integrate
interpolate
intkon1
intkon2
intkon3
intkon4
intvector
kern
mask
matlab_c
postproc
pot
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

