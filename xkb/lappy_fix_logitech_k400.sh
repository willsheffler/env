ids = xinput list | sed -n "s/.*Logitech\ K400.*id=\([0-9]*\).*pointer.*/\1/p"
for id in $ids; do 
    file = $HOME/env/lappy_logitech_k400.xkb
    echo "loading $file to device_id $id"
    xkbcomp -i $id $file
done
