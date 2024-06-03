#!/bin/bash
setxkbmap -option ctrl:nocaps
if [[ $HOSTNAME == "vivid" ]]; then
    echo "================ on vivid ========================="
    xkbcomp "$HOME/env/xkb/logitech_k400_RALT_to_LWIN.xkb" "$DISPLAY"
elif [[ $HOSTNAME == "lappy" ]]; then
    echo "================ on lappy ========================="
    xkbcomp "$HOME/env/xkb/lappy_logitech_k400.xkb" "$DISPLAY"
    echo 1 | sudo tee /sys/module/hid_apple/parameters/swap_opt_cmd
    # k400ids=$(xinput list | sed -n "s/.*Logitech\ K400.*id=\([0-9]*\).*pointer.*/\1/p")
    # for k400id in $k400ids; do
        # echo "k400 k400id: $k400id"
        # xkbcomp -i "$k400id" "$HOME/env/xkb/lappy_logitech_k400.xkb" "$DISPLAY"

    # done
    # internal_ids=$(xinput list | sed -n "s/.*Apple\ Internal\ Keyboard.*id=\([0-9]*\).*keyboard.*/\1/p")
    # for internal_id in $internal_ids; do
        # echo "k400 internal_id: $internal_id"
        # xkbcomp -i "$internal_id" "$HOME/env/xkb/lappy_internal_orig.xkb" "$DISPLAY"
    # done
elif [[ $HOSTNAME == "phat" ]]; then
    echo "================ on phat =========================="
    xkbcomp "$HOME/env/xkb/dell7720.xkb" "$DISPLAY"
else
    echo "don't know about hostname $HOSTNAME"
fi
