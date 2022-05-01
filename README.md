# MaxInfBRkNN
The source codes for "Maximizing the Influence of Bichromatic Reverse k Nearest Neighbors in Geo-Social Networks" 
Please refer to paper: https://arxiv.org/abs/2204.10203 for details.

## Running Environment
a 64-bit Linux-based OS;

## Files 
All of the C++ source code can be found in "MaxInfBRkNN" directory, the datasets can be downloaded from https://pan.baidu.com/s/1sEZ_2OWEFqaHsII1JB5ziQ, and please put all the files of the corresponding dataset into the subdirectory of "data" with the same dataset name;

## Setup
Open the "config.h" file in MaxInfBRkNN directory, and select the dataset (i.e., Gowalla or LasVegas) for testing;

## How to run 
- **Step 1: clean the source datasets**      
```shell
./FileName  -e pre_process     //FileName is name of complied executable file
```
- **Step 2: compute distance lower and upper bounds matrixs(requried in the basline solution)** 
```shell     
./FileName  -e disbound       
```
- **Step 3: modify the keyword frequency under Zipf distributions**
```shell      
./FileName  -e text           
```
- **Step 4: process the social network struture**      
```shell 
./FileName  -e social          
``` 
- **Step 5: construct the hierarchical structure of the Gimtree**      
```shell 
./FileName  -e gimtree         
``` 
- **Step 6: construct the NVDs in our index**   
```shell 
./FileName   -e nvd         
```

- **Step 7: test the effectiveness of POI selection**   
```shell 
./FileName   -e effective -t g  -s (appro|rele|influ|random) -g 1      
```
- **Step 8: test the performance of our proposed methods**   
```shell 
./FileName   -e performance -t (k|b|pc|alpha)  -s (base|appro|heuris) -g 1   //test the performance of our methods (base|appro|heuris) when varying paramter settings   
```

## Remark
If there are any problems, please contact jpf@zju.edu.cn
