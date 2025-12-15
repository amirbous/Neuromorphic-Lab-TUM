#!/usr/bin/env bash

Base_problem_name=$1
number_of_problems=$2
number_of_problems=$((number_of_problems - 1))
clean_directory=$3



echo "   problem_name,  n_vertices,  num_non_zeros,   max_e_length,    l2_res_norm"

for i in $(seq 0 ${number_of_problems}); do
    complete_name=${Base_problem_name}${i}
    ./poissfem ${complete_name} 1
    python3 solve_csr_scipy.py ${complete_name}
    ./poissfem ${complete_name} 0
    if [ "$clean_directory" = 1 ]; then
        rm ${complete_name}_mtx.txt ${complete_name}_rhs.txt ${complete_name}_sol.txt
    fi
    
done
