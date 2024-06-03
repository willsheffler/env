#!/bin/sh

# The remote control device (USB keyboard)
# has keys that overlap with standard keyboards.
# Therefore we map the offending keys here.

# re-execute this script as "$xuser" if that's not root,
# in case the root user is not permitted to access the X session.
xuser=$(ps -C Xorg -C X -ouser=)
[ $(id -un) = root ] && [ "$xuser" != root ] && exec su -c $0 "$xuser"
export DISPLAY=:0

echo "xuser = " $xuser

# First find the X input ID corresponding to the
# keyboard we're interested in.
kbd_ids() {
  keyboard="$1"
  xinput list |
  sed -n "s/.*$keyboard.*id=\([0-9]*\).*pointer.*/\1/p"
}

# In our case the USB:ID was shown in the `xinput list` output,
# but this is unusual and you may have to match
# on names or even correlate with /dev/input/by-id/*
k400ids=$(kbd_ids 'K400')
[ "$k400ids" ] || exit

echo "k400ids = " $k400ids

# 133 = Super_L
# 108 = ISO_Level3_Shift
# Write out the XKB config to remap just
# the keys we're interested in
mkdir -p /tmp/xkb/symbols
cat >/tmp/xkb/symbols/custom <<\EOF
xkb_symbols "logitech_k400" {
    key <RALT>   { [ LWIN ] };
};
EOF

# Amend the current keyboard map with
# the above key mappings, and apply to the particular device.
# Note xkbcomp >= 1.2.1 is needed to support this
for k400id in $k400ids; do 
    echo "doing xkb magic for $k400id"
    setxkbmap -device $k400id -print |
    sed 's/\(xkb_symbols.*\)"/\1+custom(logitech_k400)"/' |
    xkbcomp -I/tmp/xkb -i $k400id -synch - $DISPLAY 2>/dev/null
done

