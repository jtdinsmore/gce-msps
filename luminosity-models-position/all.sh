#!/bin/bash

if [ $# != 1 ]
then
    echo "You must supply a second argument"
    exit
fi

if [ "$1" == "hist" ]
then
    ./calc hist powerlaw &
    ./calc hist lognormal &
    ./calc hist ploeg &
    ./calc hist gautam &
    ./calc hist nptf &
elif [ "$1" == "count" ]
then
    ./calc count powerlaw &
    ./calc count lognormal &
    ./calc count nptf &
else
    echo "The argument $1 is neither hist nor count."
    exit
fi
