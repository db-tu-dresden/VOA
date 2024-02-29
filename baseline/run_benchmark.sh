
mem_node=0
hash_table_mib=1
probe_gib=1
thread_count=1
shift_size=9
collisions_count=0
repeats=7

while getopts h:p:m:t:s:c:r: flag
do 
    case "${flag}" in
        h) hash_table_mib=${OPTARG};;
        p) probe_gib=${OPTARG};;
        m) mem_node=${OPTARG};;
        t) thread_count=${OPTARG};;
        s) shift_size=${OPTARG};;
        c) collisions_count=${OPTARG};;
        r) repeats=${OPTARG};;
    esac
done

numactl -N 0 -m $mem_node ./a.out -kb $hash_table_mib -p $probe_gib -t $thread_count -s $shift_size -c $collisions_count -r $repeats