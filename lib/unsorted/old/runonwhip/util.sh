
function whipkillall 
{
for i in 01 02 03 04 05 06 07 08 09 10 11
do
   #echo $i $k
    pid=`ssh whip$i ps auxf | grep sheffler@notty | grep -v grep| awk '{print($2)}'`
    ssh whip$i kill $pid
    ssh whip$i killall rosetta_features.gcc
    ssh whip$i killall rosetta_features.gccdebug
#  ssh whip$i killall rosetta.gcc
#  ssh whip$i killall dssp
done
}

function whiplistproc 
{
    tot=0
    for i in 01 02 03 04 05 06 07 08 09 10 11
      do
      tmp=`ssh whip$i 'ps aux | egrep "rosetta|dssp" | grep sheffler | grep -v grep | wc -l'`
      echo $i $tmp
      tot=`expr $tot + $tmp`
    done
    echo $tot
}

