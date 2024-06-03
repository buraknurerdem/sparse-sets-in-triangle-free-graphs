# Sparse Sets in Triangle-free Graphs

This is the source code used in the study of Sparse Sets in Triangle-free Graphs. To run this code, you need to download nauty and compile it to work with thread-local storage. Check nauty documentation for further information.

### The Extremal Graphs

The extremal graphs for some defective Ramsey numbers are available in this repository. They are generated with the code in `main.cpp`. The graphs are provided in the adjacency matrix format. One can inspect the graphs visually using the graph drawing tools on the web that allows the input of an adjacency matrix. The name format of the files are explained below.

``T_k(j)_x.txt`` includes the graphs of highest order that do not include triangles nor k-sparse j-sets. There are x many of them.

``R_k(i,j)_x.txt`` includes the graphs of highest order that do not include triangles, k-dense i-sets nor k-sparse j-sets. There are x many of them. 
