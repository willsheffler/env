# all this is on whip01:/scratch

cd /dump/kira/
for i in `/bin/ls`; do
    if [ -d $i ]; then
	echo $i
    fi
done
