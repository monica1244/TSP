## Traveling Salesman Problem
### Course Project for CSE 6140: Algorithms (Prof. Xiuwei Zhang) @ Georgia Institute of Technology

### 4 Algorithms for TSP
In this repo, we provide 4 ways of solving the TSP problem by using Branch and Bound, Genetic Algorithm, Simulated Annealing and Approximation approaches.

### Folder structure
There are two main folders and the make file:
- code
	- bnb.cpp                               branch and bound code
	- bnb.h                                 branch and bound header file
	- genetic_algo.cpp                      local search 2 code
	- genetic_algo.h                        local search 2 header file
	- prim_dfs.cpp                          approximation algo code
	- prim_dfs.h                            approximation algo header file
	- simulated_annealing.cpp               local search 1 code
	- simulated_annealing.h                 local search 1 header file
	- TSP-main.cpp                          tsp-main calls functions from the 4 algos codes
- output
- main                                          make file
- DATA

The output folder containts all .sol and .trace files.

### Usage
There are comman line arguments for providing the input file name, algorithm type, time limit and random seed.<br>
Example of how to compile the make file:<br>
./main -inst Cincinnati.tsp -alg LS2 -time 13 -seed 3<br><br>

(ignore the warnings)

### Contributors
Sahith Dambekodi<br>
Monica Gupta<br>
Sanmeshkumar Udhayakumar<br>
Christopher Fleisher

This repo hosts the final version of the code hosted at: https://github.com/SND96/TSP-Project
