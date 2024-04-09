benchmark="baseline_ht_size"
repeats=7
node=0
probe=4
collisions_count=0
shift_size=9
thread_count=1
file_override=0
file_override_name=""

while getopts c:p:m:t:s:r:f: flag
do 
    case "${flag}" in
        c) collisions_count=${OPTARG};;
        p) probe=${OPTARG};;
        m) node=${OPTARG};;
        t) thread_count=${OPTARG};;
        s) shift_size=${OPTARG};;
        r) repeats=${OPTARG};;
        f) file_override=1 
           file_override_name=${OPTARG};;
    esac
done

memory="RAM"
if [ node >= 8]
then 
    memory="HBM"
fi

datetime=`date +"%Y-%m-%d %T"`
datetime=${datetime//[\ ]/_}
datetime=${datetime//[:]/-}

filename="$benchmark"_"$datetime"

if [ file_override > 0 ]
then
    filename="$file_override_name"
fi

filename="$filename"_"$memory"

./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 0.25;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 0.5;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 1;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 2;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 22;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 105;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 128;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 512;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 1024;