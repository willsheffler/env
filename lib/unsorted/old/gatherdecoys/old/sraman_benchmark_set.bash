#!/bin/bash

prots='1abo 1c9o'

prots='1abo  1c9o  1dtd  1ig5  1lmb  1ptq  1ubi 
       1ail  1cc8  1dxg  1iib  1lou  1r69  1utg
       1aiu  1cei  1elw  1isu  1lz1  1scj  1vcc
       1b3a  1csp  1erv  1kpe  1mzm  1tif  1vif
       1bkr  1ctf  1ew4  1kte  1orc  1tig  1wap
       1bq9  1dan  1gvp  1kve  1pcf  1tuc  1who
       1c8c  1dhn  1hz6  1lis  1pgx  1tul'

for i in $prots; do
    ls $i
done

for i in $prots; do
  tar -xf $i.tar
  if [ -d farlx ]; 
      then
      echo $i `/bin/ls farlx/ | wc -l`
      mv farlx $i
  else
      echo $i farlx missing
  fi
done

for i in $prots; do
    echo $i `'ls' $i | wc -l`
done
