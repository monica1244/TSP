## Traveling Salesman Problem
#### Course Project for CSE 6140: Algorithms (Prof. Xiuwei Zhang) @ Georgia Institute of Technology

### 4 Algorithms for TSP
In this repo, we provide 4 ways of solving the TSP problem by using Branch and Bound, Genetic Algorithm, Simulated Annealing and Approximation approaches.

### Folder structure
There are two main folders and the make file:
- code<br>
	- bnb.cpp:                               branch and bound code<br>
	- bnb.h:                                branch and bound header file<br>
	- genetic_algo.cpp:                      local search 2 code<br>
	- genetic_algo.h:                        local search 2 header file<br>
	- prim_dfs.cpp:                          approximation algo code<br>
	- prim_dfs.h:                            approximation algo header file<br>
	- simulated_annealing.cpp:               local search 1 code<br>
	- simulated_annealing.h:                 local search 1 header file<br>
	- TSP-main.cpp:                          tsp-main calls functions from the 4 algos codes<br>
- output<br>
- main                                          make file<br>
- DATA<br>

The output folder containts all .sol and .trace files.

### Usage
There are command line arguments for providing the input file name, algorithm type, time limit (optional) and random seed (optional).<br><br>
Example of how to compile the make file:<br>
./main -inst Cincinnati.tsp -alg LS2 -time 13 -seed 3<br>

(ignore the warnings)

### Contributors
Sahith Dambekodi @ https://github.com/SND96<br>
Monica Gupta @ https://github.com/monica1244<br>
Sanmeshkumar Udhayakumar @ https://github.com/sanmesh1<br>
Christopher Fleisher @ https://github.com/gravaman

This repo hosts the final version of the code hosted at: https://github.com/SND96/TSP-Project
