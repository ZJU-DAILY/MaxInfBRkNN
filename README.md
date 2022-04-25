# MaxInfBRkNN
The source codes for "Maximizing the Influence of Bichromatic Reverse k Nearest Neighbors in Geo-Social Networks" 
Please refer to paper: https://arxiv.org/abs/2204.10203 for details.

Running Environment: a 64-bit Linux-based OS;

Files: All of the C++ source code can be found in the top-level MaxInfBRkNN directory, the source (cleaned) datasets can be found in the top-level data directory;

Setup: Open the "config.h" file in MaxInfBRkNN directory, and select the dataset (i.e., Gowalla or LasVegas) for testing;

How to run (FileName is name of complied executable file): 
1.   ./FileName  -e pre_process   # Clean the source datasets 
2.   ./FileName  -e disbound      # Compute distance lower and upper bounds(which optimizes our basline solution) 
3.   ./FileName  -e text          # modify the keyword frequency under Zipf distributions
4.   ./FileName  -e text          # process the social network struture
5.   ./FileName  -e gimtree       # construct the major structure of the Gimtree
6.   ./FileName  -e nvd           # construct the NVDs in our index
7.   ./FileName  -e effective -t g  -s (appro|rele|influ|random) -g 1       # test the effectiveness of POI selection policies
8.   ./FileName  -e performance -t (k|b|pc|alpha)  -s (base|appro|heuris) -g 1  # test the performance of our methods (base|appro|heuris) when varying paramter settings

Remark: If there are any problems, please contact jpf@zju.edu.cn
