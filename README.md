# bitonic_sort
Batcher bitonic mergesort.

https://en.wikipedia.org/wiki/Bitonic_sorter

How to run and debug the program, read in the file header.

Реализована параллельная битонная сортировка в условиях ограниченной памяти.

Общая используемая память: 2 длины входного массива.

На каждом процессе не больше 2 длин входного массива деленного на количество процессов.

На одном процессе весь массив не собирается.

Правильность работы можно проверить поставив PRINT_RESULTS 1.
