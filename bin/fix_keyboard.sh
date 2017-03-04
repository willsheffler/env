#!/bin/bash
setxkbmap -option ctrl:nocaps
xkbcomp $HOME/env/logitech_k400_RALT_to_LWIN.xkb $DISPLAY
