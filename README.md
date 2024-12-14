# Accelerating Shortest Path Counting on Road Networks
This repository contains the algorithm implementations described in the paper submitted to ICDE 2025 Round 2:  *Accelerating Shortest Path Counting on Road Networks*.

# Run code
To run **CTL-Construct** and **CTL-Query**, enable (do not comment out) line 9 in the **road_network.h** file.

To run **CTLS\*-Construct**, comment out line 9 and enable line 11 in the **road_network.h** file.


`make`



### 1. Index construction 

./index graph_dir index_name

`./index dataset/d.NY.gr idx_NY`

### 2. Query 

./query index_name output_format num_of_files

`./query idx_NY dataset/NY_query dataset/ans.NY 2`

# dataset
We use distance graph datasets available from the following source:  
[https://www.diag.uniroma1.it/challenge9/download.shtml](https://www.diag.uniroma1.it/challenge9/download.shtml)


# Reference
code extend from: [road-networks](https://github.com/henningkoehlernz/road-networks)