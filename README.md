# bitonic_sort
Batcher bitonic mergesort.
Language: C.

https://en.wikipedia.org/wiki/Bitonic_sorter

To run and debug the program, check the file (.c) header.

Implemented parallel bitonic sorting under limited memory conditions.

Total memory used: 2 lengths of the input array.

Each process has at most 2 lengths of the input array divided by the number of processes.

On one process the entire array is not collected.

Correction of the calculations can be checked by setting PRINT_RESULTS 1.
