#!/bin/bash

alias rm='rm -i'
alias ls='ls -B --color'
alias more=less
alias ll='ls -B --color -l'
alias la='ls -B --color -a'
alias lla='ls -B --color -a -l'
alias hn=hostname
alias l=less

alias enscript2='enscript -fCourier7 --pretty-print --color -2r -DDuplex:true -DTumble:false --margins=30:30:30:30'
alias en2f='enscript2 --pretty-print=fortran'
alias en2p='enscript2 --pretty-print=python'
alias en2c='enscript2 --pretty-print=cpp'

alias q="condor_q | grep sheffler"
