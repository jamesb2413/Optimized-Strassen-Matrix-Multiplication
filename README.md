A project for Harvard's Data Structures and Algorithms course. Simulating Kruskal's Minimum Spanning Tree algorithm on random, 
complete, undirected graphs. A graph with n vertices is complete if all n choose 2 pairs of vertices are edges in the graph.

# Overview
Strassen’s recursive divide and conquer matrix multiplication algorithm for n by n matrices is asymptotically faster than the 
conventional 
$O(n^3)$ algorithm. Although Strassen's is faster for sufficiently large values of n, the conventional algorithm is faster 
for small values of n.

![Alt text](https://www.geeksforgeeks.org/wp-content/uploads/strassen_new.png)

At some point in Strassen's recursion, at some point in the recursion, once the matrices are small enough, we want to switch from recursively 
calling Strassen’s algorithm and just do a conventional matrix multiplication because conventional matrix multiplication is faster up to some reasonable size.
That is, instead of recursing all the way down to a base case of a 1 by 1 matrix, we will determine a larger matrix size base case where we switch
from Strassen's to the conventional algorithm. This "cross-over point" is defined here as the value of n below which we switch to the conventional
algorithm. 

The goal of this project is to experimentally determine the optimal cross-over point and use that to implement an optimal matrix multiplication
algorithm. Considering memory management was also important in writing an efficient algorithm.

# Code Layout
In Matrix.java, a square matrix data type is defined, taking a dimension input on construction. The `randLoad()` function fills a matrix with 
random integer values in a given range. 

In Strassen.java are implementations of the `standard()` and `strassen()` matrix multiplication algorithms.

I optimized memory management in Strassen’s by avoiding excessive memory allocation and deallocation and by avoiding copying large blocks of data unnecessarily. 
Instead of initializing 8 n/2 matrices A,...,H and then adding and subtracting them into the sum matrices, I only created 14 matrices and recycled them. Referencing the lecture 10 notes for Strassen’s algorithm, I first stored A, F − H, A + B, H, C + D, etc. in the 14 matrices, then overwrote 7 of the matrices to become P1, P2, ..., P7. Then, I again over- wrote 4 matrices to become AE + BG, AF + BH, etc. By only storing 14 matrices and overwriting matrices that were no longer needed, I reduced the amount of memory allocation and copying. Luck- ily, each step of the algorithm required sequential operations on some matrices that were not used in future steps, so this optimization was possible. Since regular Strassen’s was shown in Lecture 10 to be Θ(nlog 7), which is o(n3), our implementation is O(nlog 7) since we improve upon the Strassen’s demonstrated in Lecture using our cross-over points.


# Experiments and Discussion
To find the cross-over point experimentally, I randomly chose seven Matrix sizes between 1000 and 2000, then recorded the time to run 
Strassen’s with about ten different cross-over point values ($n_0$) (between 50 and 600) for each matrix size. For each $n_0$ value, I 
ran Strassen’s 5 times and took the average run time. I only considered sizes between 1000 and 2000 because I observed that smaller matrix 
multiplications did not have significant differences in run time between different $n_0$ values. I loaded the matrices with randomly generated zeroes and ones
for my experiments. I also tested larger randomly chosen positive integers, up to 100, and observed that there was no 
difference in run times compared to simply using 0 and 1. 

For large matrix sizes, the differences in speed resulting from different cross-over points were obvious. Often, there was a ”sweet spot,” a
range of $n_0$ values that resulted in the fastest run times, such that smaller and larger $n_0$ values outside of this range had slower 
run times. This makes sense because each matrix size makes large recursive jumps, so we find a range of possible optimal $n_0$ values from 
examining a single matrix size, rather than a single optimal $n_0$. For example, multiplying matrices of size 1224 will require recursive 
multiplications of sizes 612, 306, 153, etc., but there will never be a recursive multiplication of matrices of size 306 < n < 612. 
Therefore, if we find that this matrix multiplication is fastest when $n_0 = 500$, it should be similarly fast with $n_0 = 611$ or $n_0 = 375$
because any $n_0$ in the range $306 <= n < 612$ will result in the same behavior, namely switching to conventional multiplication on size 
306 matrices. By testing many matrix sizes, I was able to narrow down the sweet spot to find a small range of optimal $n_0$ values for all matrix sizes.
All matrix multiplications contained [256, 347] in their optimal $n_0$ range. Therefore, we can comfortably 
say that the optimal experimental $n_0$ is about 260.


Some matrix sizes had much larger optimal ranges, while others had very small optimal ranges. This may be because there
is a different optimal cross-over when the last recursive matrix size is odd vs. even. I conjecture that the large ranges of optimal $n_0$
occur for matrix sizes that land on both even and odd sizes in recursion, such that an optimal odd $n_0$ is struck and then an optimal even 
$n_0$ is struck. It was difficult to experimentally parse out the different $n_0$ for odd and even matrices since our padding strategy 
leads odd matrices to recurse on even matrices; therefore, recursions on odd matrices generally show a wide range of experimentally 
optimal cross-over points.


