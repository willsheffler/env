function runwhipjob 
{
    echo running on ${whips[whip]}
    ssh -x ${whips[whip]} "source ~/scripts/gatherdecoys/run_score_energy_polar.sh $1 $2 $3" &
    tmp=`expr $whip + 1`
    whip=`expr $tmp % ${#whips[*]}`
}

function whipkillall 
{
for i in 01 02 03 04 05 06 07 08 09 10 11
do
   #echo $i $k
    pid=`ssh whip$i ps auxf | grep sheffler@notty | grep -v grep| awk '{print($2)}'`
    ssh whip$i kill $pid
    ssh whip$i killall rosetta.gcc
    #ssh whip$i killall rosetta.gccdebug
#  ssh whip$i killall rosetta.gcc
#  ssh whip$i killall dssp
done
}

function whiplistproc 
{
    tot=0
    for i in 01 02 03 04 05 06 07 08 09 10 11
      do
      tmp=`ssh whip$i 'ps aux | egrep "rosetta.*gcc*|dssp" | grep sheffler | grep -v grep | wc -l'`
      echo $i $tmp
      tot=`expr $tot + $tmp`
    done
    echo $tot
}

scp ~/scripts/gatherdecoys/run_score_energy_polar.sh whip:scripts/gatherdecoys/

whips=(whip01 whip02 whip03 whip04 whip05 whip06  whip08 whip09 whip10 whip11)
whip=0

if [ ]; then

dir=~/data/decoys
set=kira_bench
prots="1acp  1b72  1ig5  1r69  1tig  1vii 
       1aa3  1ail  1bf4  1pgx  1tif  1utg  256b" # 39  1a32
for i in $prots; do 
    runwhipjob $dir $set $i > run_whip_job_${set}_${i}.log 2> run_whip_job_${set}_${i}.err
done


dir=~/data/decoys
set=sraman_nmr
prots="1fvk  1hb6  1hoe  1i27  1who  2int  3il8  4rnt" # 65
#prots=" 2int 4rnt" # 65
for i in $prots; do  
    runwhipjob $dir $set $i > run_whip_job_${set}_${i}.log 2> run_whip_job_${set}_${i}.err
done

dir=~/data/decoys
set=kira_hom
prots="1a4p  1agi  1ar0  1b07  1bmb  1cfy  1d2n  1erv  1pva 
       1acf  1ahn  1aw2  1be9  1btn  1crb  1dol  1ig5  1qav" # 57
#prots="1acf 1ahn 1aw2 1cfy 1crb 1erv"
for i in $prots; do  
    runwhipjob $dir $set $i > run_whip_job_${set}_${i}.log 2> run_whip_job_${set}_${i}.err
done

dir=~/data/decoys
set=bqian_caspbench
prots="t196  t199  t205  t206 t223  t224  t232  t234  t249" # 25
#prots="t199  t205  t206 t223  t234" # 25
for i in $prots; do  
    runwhipjob $dir $set $i > run_whip_job_${set}_${i}.log 2> run_whip_job_${set}_${i}.err
done

fi

if [ asdf ]; then

dir=~/data/decoys
set=phil_homolog
prots="1af7  1csp  1di2  1mky  1n0u  1ogw  1shf  1tig
       1b72  1dcj  1dtj  1mla  1o2f  1r69  1tif  2reb" #16
#prots="1af7  1csp  1di2  1mky  1shf  1tig
#       1b72  1dcj  1dtj  1mla  1r69  1tif  2reb" #16
#prots="1n0u 1ogw"
for i in $prots; do  
    runwhipjob $dir $set $i > run_whip_job_${set}_${i}.log 2> run_whip_job_${set}_${i}.err
done

fi
