#!/bin/bash

perl -i -pe 's/\n//g' list.h
sed -i 's/{//g' list.h
sed -i 's/}//g' list.h
sed -i 's/, /,/g' list.h
sed -i 's/, /,/g' list.h
sed -i 's/, /,/g' list.h
sed -i 's/( /(/g' list.h
sed -i 's/( /(/g' list.h
sed -i 's/( /(/g' list.h
sed -i 's/\/ /\//g' list.h
sed -i 's/\/ /\//g' list.h
sed -i 's/\/ /\//g' list.h
sed -i 's/\[ /\[/g' list.h
sed -i 's/\[ /\[/g' list.h
sed -i 's/\[ /\[/g' list.h
sed -i 's/ S/*S/g' list.h
sed -i 's/,0,/, 0,/g' list.h
sed -i 's/\//.0\//g' list.h
sed -i 's/\].0/\]/g' list.h
sed -i 's/\]/)/g' list.h
sed -i 's/Sqrt\[/sqrt(/g' list.h
perl -i -pe 's/,/,\n/g' list.h
