benchmark="baseline_scale"  # benchmark name
repeats=7                       # test repeats
node=0                          # memory node
probe=4                         # probe amount in gib
hash_table_size=32 
collisions_count=0    
thread_count=1
file_override=0
file_override_name=""
while getopts h:c:p:m:t:r:f: flag
do 
    case "${flag}" in
        c) collisions_count=${OPTARG};;
        h) hash_table_size=${OPTARG};;
        p) probe=${OPTARG};;
        m) node=${OPTARG};;
        t) thread_count=${OPTARG};;
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

for i in {0..9}
do
    ./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $i -r $repeats -f $filename -h $hash_table_size;
done
