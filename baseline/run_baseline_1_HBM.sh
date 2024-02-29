repeats=7
node=8
probe=4
collisions_count=0
shift_size=9
thread_count=1
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -h 1;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -h 2;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -h 22;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -h 105;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -h 128;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -h 512;
./run_benchmark.sh -m $node -p $probe -t $thread_count -c $collisions_count -s $shift_size -r $repeats -h 1024;