#!/bin/bash

for i in `/bin/ls`; do 
    if [ -d $i ]; then 
	echo $i `tail -n2 $i/log/make.fasc.files.Rout`; 
    fi ; 
done
echo `tail -n1 gather.log`
