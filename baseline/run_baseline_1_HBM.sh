benchmark="baseline_1_HBM"
repeats=7
node=8
probe=4
collisions_count=0
shift_size=9
thread_count=1

datetime=`date +"%Y-%m-%d %T"`
datetime=${datetime//[\ ]/_}
datetime=${datetime//[:]/-}

filename="$benchmark"_"$datetime"
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 0.25;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 0.5;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 1;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 2;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 22;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 105;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 128;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 512;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -f $filename -h 1024;