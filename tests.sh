
for size in 600 1200 2000
do
	for _ in {0..29}
	do
		echo Executing sequential code for N=$size seed=42 
		(time ./jacobiseq $size 42 5) 2>&1 | cat
		echo 

		for threads in 4 8 16
		do
			echo Executing parallel code for N=$size and Num_threads=$threads
			(time ./jacobipar $size $threads 42 5) 2>&1 | cat
			echo
		done
	done
done
