benchmark="baseline_threads"  # benchmark name
repeats=7                       # test repeats
node=0                          # memory node
probe=4                         # probe amount in gib
hash_table_size=32 
collisions_count=0
file_override=0
file_override_name=""
max_threads=12
shift_size=9
while getopts h:c:p:m:t:r:f: flag
do 
    case "${flag}" in
        c) collisions_count=${OPTARG};;
        h) hash_table_size=${OPTARG};;
        p) probe=${OPTARG};;
        m) node=${OPTARG};;
        t) max_threads=${OPTARG};;
        s) shift_size=${OPTARG};;
        r) repeats=${OPTARG};;
        f) file_override=1 
           file_override_name=${OPTARG};;
    esac
done

memory="RAM"
if ((node >= 8))
then 
    memory="HBM"
fi

datetime=`date +"%Y-%m-%d %T"`
datetime=${datetime//[\ ]/_}
datetime=${datetime//[:]/-}

filename="$benchmark"_"$datetime"

if ((file_override > 0))
then
    filename="$file_override_name"
fi

filename="$filename"_"$memory"

for i in `seq 1 $max_threads`
do
    ./run_benchmark.sh -m $node -p $probe -t $i -c $collisions_count -s $shift_size -r $repeats -f $filename -h $hash_table_size;
done
