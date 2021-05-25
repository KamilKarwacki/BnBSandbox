#!/bin/bash

objCode=$1

if [[ $objCode == "Parabola" ]]
then
	echo "setting parabola as objective"
	cp Parabola2D.h objective.h
elif [[ $objCode == "XQube" ]]
then
	echo "setting x³ function as objective"
	cp XQube.h objective.h
elif [[ $objCode == "sincos" ]]
then
	echo "setting sin(x)*cos(y) function as objective"
	cp sincos.h objective.h
fi

