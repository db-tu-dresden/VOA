benchmark="baseline_2-1_HBM"
repeats=7
node=8
probe=4
hash_table_size=64
shift_size=9
thread_count=1

datetime=`date +"%Y-%m-%d %T"`
datetime=${datetime//[\ ]/_}
datetime=${datetime//[:]/-}

filename="$benchmark"_"$datetime"
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 1 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 2 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 3 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 4 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 5 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 6 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 7 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c 8 -s $shift_size -r $repeats -f $filename -h $hash_table_size;
