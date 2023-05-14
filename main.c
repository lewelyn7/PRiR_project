#include <stdlib.h>
#include <stdio.h>
#include "omp.h"
#include <time.h>


void present_result(int mode, int threads, long  size, long  k, int result, double time){
    FILE *filePointer ;
    filePointer = fopen("results.csv", "a") ;
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


int parallel_quickselect_no_alloc(int a[], long  low, long  index_to_select,
		long  high, const long  num_threads, bounds t_bounds[], long  p_indexes[]) {
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

	if(index_to_select >= size)
		return -1;

	if(num_threads*2>=size)
		return quickselect(a, low, index_to_select, high);

	for(i=0;i<num_threads;i++) {
		t_bounds[i].low = low + i * size / num_threads;
		t_bounds[i].high = i < num_threads - 1 ? low + (i+1) * size / num_threads - 1 : high;
	}

	while(1) {

		pivot_index_on_current_partitions = rand() % size;
		for(i=0,pivot_index=0;pivot_index_on_current_partitions >= 0; i++) 
			pivot_index_on_current_partitions -= t_bounds[i].high - t_bounds[i].low + 1;
		
		i--;
		pivot_index_on_current_partitions += t_bounds[i].high - t_bounds[i].low + 1;
		pivot_index = t_bounds[i].low + pivot_index_on_current_partitions;
		pivot = a[pivot_index];

        #pragma omp parallelfor schedule(static,1)
		for(i=0;i<num_threads;i++)
			p_indexes[i]=partition(a, t_bounds[i].low, t_bounds[i].high, pivot);


		smaller_than_pivot_count = 0;
		for(i=0;i<num_threads;i++)
			smaller_than_pivot_count+=p_indexes[i]-t_bounds[i].low;

		if(smaller_than_pivot_count == index_to_select)
			break;

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
		long  high, long  num_threads) {
	bounds *t_bounds;
	long  *p_indexes;
	int result;

	t_bounds = malloc(sizeof(bounds) * num_threads);
	p_indexes = malloc(sizeof(long ) * num_threads);
	result = parallel_quickselect_no_alloc(a, low, index_to_select,
			high, num_threads, t_bounds,p_indexes);

    free(t_bounds);
    free(p_indexes);

	return result;
}


int calc_parallel(int * arr, long  size, long  k, int threads){
    double start_time, end_time;
    start_time =  omp_get_wtime();
    int parallel_result = 0;

    parallel_result = parallel_quickselect(arr, 0, size-k, size-1, 8);

    end_time =  omp_get_wtime();
    present_result(2, threads, size, k, parallel_result, end_time - start_time);
    return parallel_result;
}

int calc_single(int * arr, long  size, long k, int threads){
    double start_time, end_time;
    start_time =  omp_get_wtime();
    long  simple_result = parallel_quickselect(arr, 0, size-k, size-1, 1);
    end_time = omp_get_wtime();
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



void copy_array(int * arr_source, int * arr_dest, long  size){
    for(long  i = 0; i < size; i++){
        arr_dest[i] = arr_source[i];
    }
}
int main(int argc, char ** argv) {
    srand(time(NULL));
	if(argc != 3){
		printf("wrong arguments");
		exit(-1);
	}
    long  n = atol(argv[1]); // size of the array
    long  k = 345; // find the kth biggest number
    int threads = atoi(argv[2]);
    int * main_arr;
    int * working_arr;
    int single_result, parallel_result, check_result;

    main_arr = allocate_array(n);
    fill_random(main_arr, n);

    omp_set_num_threads(threads);

    working_arr = allocate_array(n);
    copy_array(main_arr, working_arr, n);
    single_result = calc_single(working_arr, n, k, threads);
    free(working_arr);

    working_arr = allocate_array(n);
    copy_array(main_arr, working_arr, n);
    parallel_result = calc_parallel(working_arr, n, k, threads);
    free(working_arr);

    working_arr = allocate_array(n);
    copy_array(main_arr, working_arr, n);
    // check_result = calc_check(working_arr, n, k, threads);
    free(working_arr);


    free(main_arr);
    if (single_result != check_result || parallel_result != check_result)
    {
        return -1;
    }
    

    return 0;
}
