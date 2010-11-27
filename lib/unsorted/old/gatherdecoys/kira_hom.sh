source /users/sheffler/scripts/gatherdecoys/gatherdecoysfunc.sh

khs='/data/kira/homology_models/homolog_vall/'
kht='/users/sheffler/data/decoys/kira_hom/'
khnat='/data/kira/homology_models/'
khp='1a4p 1acf 1agi 1ahn 1aoy 1ar0 1aw2 1b07 1b34 1be9 1bmb 1btn 1c9o 1cfy 1crb 1d2n 1dol 1e8a 1erv 1ig5 1pva 1qav 1tim'
khdt='_l3/'
khdl='none'
khst='_l3/'
khsl='nstruct|output'

function gather {
for i in $khp; do
gatherdecoys $khs/$i $kht/$i $khnat/$i/${i}_minimized_native $i $khdt $khdl $khst $khsl $khs/../$i/${i}_frags
done
}

mkdir -p $kht
gather > $kht/gather.log 2> $kht/gather.err

cp -r /data/kira/homology_models//1ahn/1ahn_minimize_native $kht/1ahn/relax_native
