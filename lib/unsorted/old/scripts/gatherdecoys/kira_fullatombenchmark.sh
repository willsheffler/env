source /users/sheffler/scripts/gatherdecoys/gatherdecoysfunc.sh

ks='/scratch/sheffler/kira_fullatombenchmark/'
kt='/users/sheffler/data/decoys/kira_fullatombenchmark/'
knat='/users/sheffler/data/relax_native/kira_fullatombenchmark/'
kp='1a32 1aa3 1acp 1ail 1b72 1bf4 1ig5 1pgx 1r69 1tif 1tig 1utg 1vii 256b'
kdt='.'
kdl='rescore'
kst='.'
ksl="'rescore|kira'"

function gather {
for prot in $kp; do
    gatherdecoys $ks/$prot $kt/$prot $knat/$prot/relax_native $prot $kdt $kdl $kst $ksl $ks/fragments
done
}

mkdir -p $kt
gather > $kt/gather.log 2> $kt/gather.err
