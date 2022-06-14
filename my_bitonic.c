/* C Program for Bitonic Sort. Note that this program
   works only when size of input and number of processes are a power of 2
   
   input: number of processes (must be a power of 2) and degree
   prints: calculating time of an array with random numbers with length = 2^degree
        set PRINT_ARRAYS 1 to print initial and sorted arrays

   run: mpicc my_bitonic.c -o my_bitonic
        mpirun -n num_processes /path/to/my_bitonic degree
   
   Sorts an array in ascending order. To sort in descending order change dir to 0
   author: Komarov Maksym (Jun 2022)
    */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

#define PRINT_ARRAYS 0

double timer_start;
double timer_end;
int process_rank;
int num_processes;
long long array_size;
long long global_array_size;

int check_inputs(int argc, char * argv[])
{    
    // Set degree from cmd
    if (argc != 2) 
    {
        printf("ERROR WRONG argc = %d \n", argc);
        return -1;
    }

    int degree = atoi(argv[1]);
    global_array_size = 1;
    for (int i = 0; i < degree; i++)
        global_array_size *= 2;

    if (degree == 0 || global_array_size < num_processes)
    {
        printf("ERROR: TOO LOW DEGREE\n");
        return -2;
    }
    
    int tmp_num_processes = num_processes;
    for (int i = 1; i < num_processes; i++)
    {
        if (tmp_num_processes % 2 != 0 && tmp_num_processes != 1)
            {
                printf("ERROR: num_processes must be a power of 2\n");
                return -3;
            }

        if (tmp_num_processes == 1)
            break;

        tmp_num_processes /= 2;
    }

    return 1;
}

int min(long long num1, long long num2) 
{
    return (num1 > num2) ? num2 : num1;
}

//If global_pointer in another process data part we need send and recv data
int needSend(long long global_pointer)
{
    // part of data located in grather process OR part of data located in less process
    if ((global_pointer >= global_array_size * (process_rank + 1) / num_processes) || 
        (global_pointer < global_array_size  * (process_rank    ) / num_processes))
        return 1;
    else
        return 0;
}

//Get rank of process that contains data part including global_pointer
int sendTo(long long global_pointer)
{
    return num_processes * global_pointer / global_array_size;
}

//get the global start index of the current process array
long long getGlobalStart()
{
    return global_array_size * process_rank / num_processes;
}


/*The parameter dir indicates the sorting direction, ASCENDING
   or DESCENDING; if (a1[i] > a2[j]) agrees with the direction,
   then a1[i] and a2[j] are interchanged.*/
void compAndSwap(long long a1[], long long a2[], long long i, long long j, int dir)
{
    if (dir==(a1[i]>a2[j]))
    {
        long long tmp = a1[i];
        a1[i] = a2[j];
        a2[j] = tmp;
    }
}
 
/*It recursively sorts a bitonic sequence in ascending order,
  if dir = 1, and in descending order otherwise (means dir=0).
  The sequence to be sorted starts at index position start_index,
  the parameter cnt is the number of elements to be sorted.*/
void bitonicMerge(long long a[], long long start_index, long long cnt, int dir)
{
    if (cnt>1)
    {   
        long long k = cnt/2;
        long long to_local_pointer = getGlobalStart();

        long long arr_start = to_local_pointer;
        long long arr_end = to_local_pointer + array_size - 1;

        //--------------------------------------CALCULATION OF A PAIR OF ARRAYS TO COMPARE AND SWAP--------------------------------------
        // CALCULATION of indexes of beginnings of arrays
        // if merge local elements:: Global index of beginings: a1) start_index and a2) start_index + k
        // if merge between processes: 
        // 1) arr_start in 1 part :: a1) arr_start   a2) start_index + k + (arr_start - start_index)        (It is bias in brackets)
        // 2) arr_start in 2 part :: a1) start_index + (arr_start - (start_index + k))  a2) arr_start       (It is bias in brackets)
        long long* a1 = (a + start_index - to_local_pointer);
        long long* a2 = (a + start_index + k - to_local_pointer);
        //printf("r: %d, k = %lld a1 starts: %lld, a2 atarts: %lld \n", process_rank, k, start_index - to_local_pointer, start_index + k - to_local_pointer);
        
        // if it is not local (not 1 process) bitonicMerge:
        // 1) if the process does not include the second part of the global array to merge
        if (arr_end < start_index + k)
        {
            // 1) arr_start in 1 part :: Global index of the beginings:
            // a1) arr_start   a2) start_index + k + (arr_start - start_index)        (It is bias in brackets)
            long long a2_global_start = arr_start + k;
            a1 = a; // local start index is 0
            
            //printf("\nSend from %d to %d\n", process_rank, sendTo(a2_global_start));
            // Send array to paired process
            MPI_Send(
                a1,
                array_size,                      // enire array
                MPI_LONG_LONG,
                sendTo(a2_global_start),         // paired process from second half of tmp global array (with cnt elements)
                cnt,                             // tag
                MPI_COMM_WORLD                   // default comm.
            );

            //printf("\nRecv in %d from %d\n", process_rank, sendTo(a2_global_start));
            // Receive array from paired process
            a2 = (long long*)calloc(array_size, sizeof(long long));
            MPI_Recv(
                a2,                              // data from paired process
                array_size,                      // one data item
                MPI_LONG_LONG,               
                sendTo(a2_global_start),         // paired process from second half of tmp global array (with cnt elements)
                cnt,                             // tag
                MPI_COMM_WORLD,                  // default comm.
                MPI_STATUS_IGNORE                // ignore info about message received
            );
        }
        // if it is not local (not 1 process) bitonicMerge:        
        // 2) if the process does not include the first part of the global array to merge
        if (arr_start >= start_index + k)
        {         
            // 2) arr_start in 2 part :: Global index of the beginings:
            // a1) start_index + (arr_start - (start_index + k))  a2) arr_start       (It is bias in brackets)
            long long a1_global_start = arr_start - k;
            a2 = a; // local start index is 0

            //printf("\nRecv in %d from %d\n", process_rank, sendTo(a1_global_start));
            // Receive array from paired process
            a1 = (long long*)calloc(array_size, sizeof(long long));
            MPI_Recv(
                a1,                              // data from paired process
                array_size,                      // one data item
                MPI_LONG_LONG,               
                sendTo(a1_global_start),         // paired process from second half of tmp global array (with cnt elements)
                cnt,                             // tag
                MPI_COMM_WORLD,                  // default comm.
                MPI_STATUS_IGNORE                // ignore info about message received
            );
            
            //printf("\nSend from %d to %d\n", process_rank, sendTo(a1_global_start));
            // Send array to paired process
            MPI_Send(
                a2,
                array_size,                      // enire array
                MPI_LONG_LONG,
                sendTo(a1_global_start),         // paired process from second half of tmp global array (with cnt elements)
                cnt,                             // tag
                MPI_COMM_WORLD                   // default comm.
            );
        }

        //----------------------------------------------COMPARE AND SWAP----------------------------------------------
        //printf("r: %d, compAndSwap part: %lld , %lld\n", process_rank, start_index, start_index + k);
        long long i;
        for (i = 0; i < min(k, array_size); i++) 
        {
            // calculate local pointers
            //printf("\n rank = %d; i = %lld. a1[i], a2[i]: %lld, %lld. To local poiner: %lld\n", process_rank, i, a1[i], a2[i], to_local_pointer);
            compAndSwap(a1, a2, i, i, dir);
        }
        
        // free memory
        if (arr_end < start_index + k)
            free(a2);
            
        if (arr_start >= start_index + k)
            free(a1);
        
        //printf("r: %d, need send: %d, %d\n", process_rank, needSend(start_index), needSend(start_index+k));//  to_local_pointer -> to_local_pointer + array_size
        // if current process array is in the first part of merge array OR
        // we need to merge local elements
        if ((arr_start >= start_index && arr_end < start_index + k) ||
            (array_size > k))
            bitonicMerge(a, start_index, k, dir);

        // if current process array is in the second part of merge array OR
        // we need to merge local elements
        start_index += k;
        if ((arr_start >= start_index && arr_end < start_index + k) ||
            (array_size > k))
            bitonicMerge(a, start_index, k, dir);
    }
}

/* This function first produces a bitonic sequence by recursively
    sorting its two halves in opposite sorting orders, and then
    calls bitonicMerge to make them in the same order */
/*
bitonicSort(entire array, define order)
{
    bitonicSort(first part, ascending order);
    bitonicSort(second part, descending order);
    bitonicMerge(entire array, define order);
}
*/
void bitonicSort(long long a[], long long start_index, long long cnt, int dir)
{
    if (cnt>1)
    {
        long long k = cnt/2;

        // if current process array part is in range(start_index, start_index + k) OR
        // we need to Sort local elements
        long long arr_start = getGlobalStart();
        long long arr_end = getGlobalStart() + array_size - 1;
        if ((arr_start >= start_index && arr_end < start_index + k) ||
            (array_size > k))
        {
            // sort in ascending order since dir here is 1
            //printf("r: %d, 1) bitonicSort: start_index: %lld, k: %lld\n", process_rank, start_index, k);
            bitonicSort(a, start_index, k, 1);
        }

        // if current process array part is in range(start_index + k, start_index + cnt) OR
        // we need to Sort local elements
        if ((arr_start >= start_index + k && arr_end < start_index + cnt) ||
            (array_size > k))
        {
            // sort in descending order since dir here is 0
            //printf("r: %d, 2) bitonicSort: start_index: %lld, k: %lld\n", process_rank, start_index+k, k);
            bitonicSort(a, start_index+k, k, 0);
        }
 
        // Will merge whole sequence in ascending order
        // since dir=1.
        //printf("r: %d, 3) bitonicMerge: start_index: %lld, cnt: %lld\n", process_rank, start_index, cnt);
        MPI_Barrier(MPI_COMM_WORLD);
        bitonicMerge(a, start_index, cnt, dir);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
 
/* Caller of bitonicSort for sorting the entire array of
   length N in ASCENDING order */
void sort(long long a[], long long N, int dir)
{
    bitonicSort(a, 0, N, dir);
}
 

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    assert(check_inputs(argc, argv) == 1);

    // Set parameters
    int degree = atoi(argv[1]);
    int dir = 1; // means sort in ascending order

    global_array_size = 1;
    for (int i = 0; i < degree; i++)
        global_array_size *= 2;

    // Initialize Array for Storing Random Numbers
    array_size = global_array_size / num_processes;
    long long* a = (long long *)malloc(array_size * sizeof(long long));
    
    // Set rand seed
    srand(process_rank + 312);
    // Complete the array by random numbers
    for (long long i = 0; i < array_size; i++) 
    {
        a[i] = rand() % global_array_size;
    }
    if (PRINT_ARRAYS)
    {
        if ((process_rank == 0) || MPI_Recv(&dir, 1, MPI_INT, process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) == MPI_SUCCESS)
        {
            fflush(stdout);
            printf("r: %d, ARRAY: ", process_rank);
            for (long long i = 0; i < array_size; i++) 
                printf("%lld ", a[i]);
            printf("\n");
            if (process_rank != num_processes - 1)
                MPI_Send(&dir, 1, MPI_INT, process_rank + 1, 0, MPI_COMM_WORLD);
        }
        printf("\n");
    }
 
    // Store start time to get sort time
    double start_time;
    MPI_Barrier(MPI_COMM_WORLD);
    if (process_rank == 0)
        start_time = MPI_Wtime();

    // dir = 1 // means sort in ascending order
    sort(a, global_array_size, dir);
    MPI_Barrier(MPI_COMM_WORLD);

    double finish_time = MPI_Wtime();
    if(process_rank == 0)
        printf("\nCalculating time: %f \n", finish_time - start_time);

    // Recv time to get calculating time

    // print arrays
    if (PRINT_ARRAYS)
    {
        if ((process_rank == 0) || MPI_Recv(&dir, 1, MPI_INT, process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) == MPI_SUCCESS)
        {
            fflush(stdout);
            printf("r: %d, Sorted array: \n", process_rank);
            for (long long i = 0; i < array_size; i++)
                printf("%lld ", a[i]);
            printf("\n");
            if (process_rank != num_processes - 1)
                MPI_Send(&dir, 1, MPI_INT, process_rank + 1, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    free(a);
    MPI_Finalize();
    return 0;
}