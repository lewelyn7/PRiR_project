#include <stdlib.h>
#include <stdio.h>
#include "omp.h"
#include <time.h>
#include <limits.h>
#include <mpi.h>


void present_result(int mode, int threads, long  size, long  k, int result, double time){
    FILE *filePointer ;
    filePointer = fopen("results_mpi.csv", "a") ;
    if ( filePointer == NULL )
    {
        printf( "results file failed to open");
    }
    else
    {
         
        fprintf(filePointer, "%d;%d;%ld;%ld;%d;%f\r\n", mode, threads, size, k, result, time);
        fclose(filePointer) ;
    }
    printf("Mode: %d threads: %d size: %ld %ldth biggest number is:  %d took: %f\r\n", mode, threads, size, k, result, time);
}
long  simplepartition(int a[], long  low, long  high) {
	int tmp;

	int pivot = a[low + rand() % (high - low + 1)];

	while (low < high) {
		while (a[low] < pivot)
			low++;

		while (a[high] > pivot)
			high--;

		if (a[low] == a[high])
			low++;
		else if (low < high) {
			tmp = a[low];
			a[low] = a[high];
			a[high] = tmp;
		}
	}

	return high;
}

int quickselect(int a[], long  low, long  index_to_select,
		long  high) {
	long  pivot_position, smaller_than_pivot_count;

	if ( low == high )
		return a[low];
	pivot_position = simplepartition(a, low, high);

	smaller_than_pivot_count = pivot_position - low;

	if (smaller_than_pivot_count == index_to_select)
		return a[pivot_position];
	else if (index_to_select < smaller_than_pivot_count)
		return quickselect(a, low, index_to_select, pivot_position - 1);
	else
		return quickselect(a, pivot_position, index_to_select - smaller_than_pivot_count, high);
}


long  partition(int a[], long  low,  long  high, long  pivot) {
	int tmp;
	const int init_high=high;
	while (low < high) {
		while (a[low] < pivot && low < high)
			low++;

		if(low==high)
			break;

		while (a[high] > pivot && high > low)
			high--;

		if(low==high)
			break;

		if (a[low] == a[high])
			low++;
		else if (low < high) {
			tmp = a[low];
			a[low] = a[high];
			a[high] = tmp;
		}
	}

	if(high==init_high && a[high] < pivot) {
		high++;
	}

	return high;
}

typedef struct bounds {
	long  low;
	long  high;
} bounds;
void print_arr(int * a, int size){
    // printf("arr: ");
    // for(int i = 0; i < size; i++){
    //     printf("%d, ", a[i]);
        
    // }
    // printf("\r\n");
}

int parallel_quickselect_no_alloc(int a[], long  low, long  index_to_select,
		long  high, const long  num_threads, bounds t_bounds[], long  p_indexes[], int * thread_buffer, int thread_buffer_size) {
	long  pivot_position, smaller_than_pivot_count;
	long  pivot;
	long  i;
	long  j,j_size;
	long  k;
	long  size = high-low+1;
	long  pivot_index_on_current_partitions;
	long  pivot_index;
	long  accum;
	long  cur_index;
	long  t_max_size, max_size;
	long  non_zero_count;
	int all_partitions_size_one_fallback[num_threads];
	const long  starting_low = low;
    MPI_Status mpi_status;

	if(index_to_select >= size)
		return -1;

	if(num_threads*2>=size)
		return quickselect(a, low, index_to_select, high);

	for(i=0;i<num_threads;i++) {
		t_bounds[i].low = low + i * size / num_threads;
		t_bounds[i].high = i < num_threads - 1 ? low + (i+1) * size / num_threads - 1 : high;
	}
    int mpi_source;
	while(1) {

		pivot_index_on_current_partitions = rand() % size;
		for(i=0,pivot_index=0;pivot_index_on_current_partitions >= 0; i++) 
			pivot_index_on_current_partitions -= t_bounds[i].high - t_bounds[i].low + 1;
		
		i--;
		pivot_index_on_current_partitions += t_bounds[i].high - t_bounds[i].low + 1;
		pivot_index = t_bounds[i].low + pivot_index_on_current_partitions;
		pivot = a[pivot_index];

		for(int i=1;i<num_threads;i++){
            int bound_size = t_bounds[i].high - t_bounds[i].low + 1;
            MPI_Send(&pivot, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
            MPI_Send(a+t_bounds[i].low, bound_size, MPI_INT, i, 1, MPI_COMM_WORLD);
        }
        p_indexes[0]=partition(a, t_bounds[0].low, t_bounds[0].high, pivot);
		for(int i=1;i<num_threads;i++){
            long rec_p_index;
            MPI_Recv(&rec_p_index, 1, MPI_LONG, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &mpi_status);
            mpi_source = mpi_status.MPI_SOURCE;
            rec_p_index += t_bounds[mpi_source].low;
            p_indexes[mpi_source] = rec_p_index;

            int bound_size = t_bounds[mpi_source].high - t_bounds[mpi_source].low + 1;
            MPI_Recv(thread_buffer, thread_buffer_size, MPI_INT, mpi_source, 3, MPI_COMM_WORLD, &mpi_status);
            // printf("0: bounds: %d %d \r\n", t_bounds[0].low, t_bounds[0].high);
            for(long j = 0; j < bound_size; j++){
                a[j+t_bounds[mpi_source].low] = thread_buffer[j];
            }
            print_arr(a, size);
        }


		smaller_than_pivot_count = 0;
		for(i=0;i<num_threads;i++)
			smaller_than_pivot_count+=p_indexes[i]-t_bounds[i].low;

		if(smaller_than_pivot_count == index_to_select){
            long finalizer = LONG_MAX;
            for(int i = 1; i < num_threads; i++){
                MPI_Send(&finalizer, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
            }
			break;
        }

		if(index_to_select < smaller_than_pivot_count) {
			for(i=0;i<num_threads;i++) {
				if(p_indexes[i] != 0)
					t_bounds[i].high = p_indexes[i] - 1;
				else
					t_bounds[i].low = t_bounds[i].high+1;
			}

		} else {
			for(i=0;i<num_threads;i++)
				t_bounds[i].low = p_indexes[i];
		}

		size = 0;
		for(i=0;i<num_threads;i++) 
			if(t_bounds[i].high >= t_bounds[i].low)
				size += t_bounds[i].high - t_bounds[i].low + 1;
		

		if(index_to_select > smaller_than_pivot_count)
			index_to_select -= smaller_than_pivot_count;

		if(num_threads>=size) {
			for(j=0,k=0; j < num_threads; j++) {
				for(i= t_bounds[j].low; i <= t_bounds[j].high; i++) {
					all_partitions_size_one_fallback[k] = a[i];
					k++;
				}
			}
            long finalizer = LONG_MAX;
            for(int i = 1; i < num_threads; i++){
                MPI_Send(&finalizer, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
            }
			return quickselect(all_partitions_size_one_fallback, 0, index_to_select, size-1);
		}

		for(i=0;i<num_threads;i++)
			if(t_bounds[i].high < t_bounds[i].low) {

				max_size=0;
				t_max_size=i;
				non_zero_count=0;
				for(j=0; j < num_threads; j++) {
					if(t_bounds[j].high > t_bounds[j].low) {
						j_size = t_bounds[j].high - t_bounds[j].low + 1;
						if(max_size<j_size) {
							max_size = j_size;
							t_max_size = j;
						}
					}
				}

				j=t_max_size;
				if(i<j) {
					t_bounds[i].low = t_bounds[j].low;
					t_bounds[i].high = (t_bounds[j].high + t_bounds[j].low) / 2;
					t_bounds[j].low = t_bounds[i].high+1;
				} else {
					t_bounds[i].high = t_bounds[j].high;
					t_bounds[j].high = (t_bounds[j].high + t_bounds[j].low) / 2;
					t_bounds[i].low = t_bounds[j].high+1;
				}
			}
	}

	return pivot;
}

int* allocate_array(long  n) {
    int* arr = (int*) malloc(n * sizeof(int));
    return arr;
}
void fill_random(int * arr, long  n){
    for(long  i = 0; i < n; i++){
        arr[i] = rand();
    }
}
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}
int find_k_biggest_serial(int * arr, long  size, long  k){
    qsort(arr, size, sizeof(int), cmpfunc);
    return arr[size-k];
}
int parallel_quickselect(int a[], long  low, long  index_to_select,
		long  high, long  num_threads, int * thread_buffer, int thread_buffer_size) {
	bounds *t_bounds;
	long  *p_indexes;
	int result;

	t_bounds = malloc(sizeof(bounds) * num_threads);
	p_indexes = malloc(sizeof(long ) * num_threads);
	result = parallel_quickselect_no_alloc(a, low, index_to_select,
			high, num_threads, t_bounds,p_indexes, thread_buffer, thread_buffer_size);

    free(t_bounds);
    free(p_indexes);

	return result;
}



int calc_parallel(int * arr, long  size, long  k, int threads, int * thread_buffer, int thread_buffer_size){
    double start_time, end_time;
    start_time =  MPI_Wtime();
    int parallel_result = 0;

    parallel_result = parallel_quickselect(arr, 0, size-k, size-1, threads, thread_buffer, thread_buffer_size);

    end_time =  MPI_Wtime();
    present_result(2, threads, size, k, parallel_result, end_time - start_time);
    return parallel_result;
}

int calc_single(int * arr, long  size, long k, int threads, int * thread_buffer, int thread_buffer_size){
    double start_time, end_time;
    start_time =  MPI_Wtime();
    long  simple_result = parallel_quickselect(arr, 0, size-k, size-1, 1, thread_buffer, thread_buffer_size);
    end_time = MPI_Wtime();
    present_result(1, threads, size, k, simple_result, end_time - start_time);
    return simple_result;
}

int calc_check(int * arr, long  size, long  k, int threads){
    double start_time, end_time;
    start_time =  omp_get_wtime();
    long  check_result = find_k_biggest_serial(arr, size, k);
    end_time =  omp_get_wtime();
    present_result(0, threads, size, k, check_result, end_time - start_time);
    return check_result;
}

void worker(int rank, int * thread_buffer, int thread_buffer_size){
    MPI_Status mpi_status; 
    int count;
    long pivot;
    while(1){
        MPI_Recv(&pivot, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &mpi_status);
        if(pivot == LONG_MAX){
            // printf("%d: exiting\r\n", rank);
            return;
        }
        // printf("%d: recieved pivot: %ld\r\n", rank, pivot); 
        MPI_Recv(thread_buffer, thread_buffer_size, MPI_INT, 0, 1, MPI_COMM_WORLD, &mpi_status);
        MPI_Get_count(&mpi_status, MPI_INT, &count);
        // printf("%d: recieved elements: %d\r\n", rank, count); 
        
        long p_index = partition(thread_buffer, 0, count-1, pivot);
        MPI_Send(&p_index, 1, MPI_LONG, 0, 2, MPI_COMM_WORLD);
        MPI_Send(thread_buffer, count, MPI_INT, 0, 3, MPI_COMM_WORLD);
    }

}

void copy_array(int * arr_source, int * arr_dest, long  size){
    for(long  i = 0; i < size; i++){
        arr_dest[i] = arr_source[i];
    }
}
int main(int argc, char ** argv) {
    srand(time(NULL));
	if(argc != 4){
		printf("wrong arguments");
		exit(-1);
	}
    long  n = atol(argv[1]); // size of the array
    int threads = atoi(argv[2]);
    long  k = atoi(argv[3]); // find the kth biggest number
	if (n % threads != 0)
	{
		n -= n % threads;
	}
    int * main_arr;
    int * working_arr;
    int single_result, parallel_result, check_result;

	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    omp_set_num_threads(threads);

	int thread_buffer_size = n / threads;
	int *thread_buffer = allocate_array(thread_buffer_size);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        main_arr = allocate_array(n);
        fill_random(main_arr, n);
    }
    if(rank == 0){
        working_arr = allocate_array(n);
        copy_array(main_arr, working_arr, n);
        single_result = calc_single(working_arr, n, k, threads, thread_buffer, thread_buffer_size);
        free(working_arr);
    }
    if(rank == 0){
        working_arr = allocate_array(n);
        copy_array(main_arr, working_arr, n);
        parallel_result = calc_parallel(working_arr, n, k, threads, thread_buffer, thread_buffer_size);
        free(working_arr);
    }else{
        worker(rank, thread_buffer, thread_buffer_size);
    }

    // working_arr = allocate_array(n);
    // copy_array(main_arr, working_arr, n);
    // // check_result = calc_check(working_arr, n, k, threads);
    // free(working_arr);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    if(rank == 0){
    free(main_arr);

    if (single_result != check_result || parallel_result != check_result)
    {
        // return -1;
    }
    }
    

    return 0;
}