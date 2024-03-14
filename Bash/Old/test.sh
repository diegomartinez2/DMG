#!/bin/bash

echo -n "Enter a number: "
read VAR

if [[ $VAR -gt 10 ]]
then
  echo "The variable is greater than 10."
fi

echo -n "Enter True: "
read myvar
myvar=`echo "True"`
case $myvar in
  (True)    echo "is true"; echo "Ture";echo "True"$myvar;;
  (False)   echo "is false";;
esac

