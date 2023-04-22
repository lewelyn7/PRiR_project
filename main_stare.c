#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

void present_result(int mode, int size, int k, int result, double time){
    printf("Mode: %d size: %d %dth biggest number is:  %d took: %f\r\n", mode, size, k, result, time);
}

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}
// https://github.com/abignoli/parallel-quickselect/blob/master/src/parqselect.hpp
int partition(int a[], unsigned long low, unsigned long high, int pivot) {
	int tmp;
	const int init_high=high;

	/*
	 * 4 Cases:
	 *
	 * 1 - pivot smaller than any element in the array: return low
	 * 2 - pivot higher than any element in the array: return high+1
	 * 3 - pivot value in the array: return rightmost occurrence of pivot value
	 * 4 - pivot value not in the array, but between max and min values in the array: return i | a[i] > pivot and a[j] < pivot for each j < i
	 */

	while (low < high) {
		while (a[low] < pivot && low < high)
			low++;

		if(low==high)
			break;

		while (a[high] > pivot && high > low)
			high--;

		if(low==high)
			break;

		if (a[low] == a[high]) // Both are equal to pivot
			low++;
		else if (low < high) {
			tmp = a[low];
			a[low] = a[high];
			a[high] = tmp;
		}
	}

	// Address Case 2
	if(high==init_high && a[high] < pivot) {
		high++;
		// high_moved = true;
	}

	return high;
}

struct bounds{
    int down;
    int up;
};
int* allocate_array(int n) {
    int* arr = (int*) malloc(n * sizeof(int));
    return arr;
}
int quickselect(int arr[], int size, int low, int high, int k, int threads) {
    int * p_indexes = allocate_array(threads);
    int * all_partitions_size_one_fallback = allocate_array(threads);
    struct bounds * bounds_arr = malloc(threads * sizeof(struct bounds));
    for(int i=0; i<threads; i++){
		bounds_arr[i].down = low + i * size / threads;
		bounds_arr[i].up = i < threads - 1 ? low + (i+1) * size / threads - 1 : high;
    }

        int pivot;
    while(1){ 
        pivot = arr[rand()%high];
        for(int i=0; i<threads; i++){
            p_indexes[i] = partition(arr, bounds_arr[i].down, bounds_arr[i].up, pivot);
        }

            // Integrate results to discover the number of elements smaller than chosen pivot
        int smaller_than_pivot_count = 0;
        for(int i=0;i<threads;i++)
            /*
                * a[p_indexes[i]] == pivot accounts for case 1,3
                * p_indexes[i] <= t_bounds[i].high accounts for case 2
                * if the conditions hold we are in case 4 and it's safe to add 1 to the size
                *
                * case 1: 0
                * case 2: size
                */
            smaller_than_pivot_count+=p_indexes[i]-bounds_arr[i].down;

        if(smaller_than_pivot_count == k)
            break;
		if(k < smaller_than_pivot_count) {
			// Pivot is in the left portions

			// Update lower and upper partition bounds for each thread
			for(int i=0;i<threads;i++) {
				/*
				 * Case 1: If pivot was smaller than any element in the partition high = low - 1 now
				 * Case 2: If pivot was higher than any element in the partition high is unchanged
				 */
				if(p_indexes[i] != 0)
					bounds_arr[i].up = p_indexes[i] - 1;
				else
/*					 If p_indexes is 0, then the partition will be surely zero-sized (since low >= 0)
					 Since high is unsigned, we cannot put -1 in it. This is relevant when low=0.
					 High should be low - 1, but high can't be -1.
					 This workaround is the equivalent as setting the partition to size 0,
					 achieving the desired result for size computation and zero-size partition fix*/
					bounds_arr[i].down = bounds_arr[i].up+1;
			}

		} else {
			// Pivot is in the right portions

			/*
			 * Case 1: If pivot was smaller than any element in the partition low is unchanged
			 * Case 2: If pivot was higher than any element in the partition low = high + 1 now
			 */
			 
			// Update lower and upper partition bounds for each thread
			for(int i=0;i<threads;i++)
				bounds_arr[i].down = p_indexes[i];
		}

		// Compute size
		size = 0;
		for(int i=0;i<threads;i++) 
			if(bounds_arr[i].up >= bounds_arr[i].down)
				size += bounds_arr[i].up - bounds_arr[i].down + 1;
		

		if(k > smaller_than_pivot_count)
			k -= smaller_than_pivot_count;

		// if size is too small fallback to standard quickselect
		if(threads>=size) {
			// Collect values in fallback array
			for(int j=0,k=0; j < threads; j++) {
				for(int i= bounds_arr[j].down; i <= bounds_arr[j].up; i++) {
					all_partitions_size_one_fallback[k] = arr[i];
					k++;
				}
			}
			return quickselect(all_partitions_size_one_fallback, threads, 0, threads-1, k, 1);
		}

		// Fix partitions of size zero
		for(int i=0;i<threads;i++)
			if(bounds_arr[i].up < bounds_arr[i].down) {
				// Thread i is assigned a partition of size 0

				// Look for thread t_max_size with partition of maximum size
				int max_size=0;
				int t_max_size=i;
				int non_zero_count=0;
                int j;
                int j_size;
				for(j=0; j < threads; j++) {
					// Finding biggest partition j to split with i
					if(bounds_arr[j].up > bounds_arr[j].down) {
						j_size = bounds_arr[j].up - bounds_arr[j].down + 1;
						if(max_size<j_size) {
							max_size = j_size;
							t_max_size = j;
						}
					}
				}

				// Note: since size > num_threads and thread i has size 0 we can guarantee that there is at least one thread with size > 1, so t_max_size > 1

				j=t_max_size;
				if(i<j) {
					// i gets left half from j
					bounds_arr[i].down = bounds_arr[j].down;
					bounds_arr[i].up = (bounds_arr[j].up + bounds_arr[j].down) / 2;
					bounds_arr[j].down = bounds_arr[i].up+1;
				} else {
					// i gets right half from j
					bounds_arr[i].up = bounds_arr[j].up;
					bounds_arr[j].up = (bounds_arr[j].up + bounds_arr[j].down) / 2;
					bounds_arr[i].down = bounds_arr[j].up+1;
				}
			}

        }
        return pivot;
    }



void fill_random(int * arr, int n){
    for(int i = 0; i < n; i++){
        arr[i] = rand();
    }
}


int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}
int find_k_biggest_serial(int * arr, int size, int k){
    qsort(arr, size, sizeof(int), cmpfunc);
    return arr[size-k];
}

int calc_parallel(int * arr, int size, int k){
    double start_time, end_time;
    start_time =  omp_get_wtime();
    int parallel_result = 0;

    parallel_result = quickselect(arr, size, 0, size-1, k, 4);

    end_time =  omp_get_wtime();
    present_result(2, size, k, parallel_result, end_time - start_time);
    return parallel_result;
}

int calc_single(int * arr, int size, int k){
    double start_time, end_time;
    start_time =  omp_get_wtime();
    int simple_result = quickselect(arr, size, 0, size-1, k, 1);
    end_time = omp_get_wtime();
    present_result(1, size, k, simple_result, end_time - start_time);
    return simple_result;
}

int calc_check(int * arr, int size, int k){
    double start_time, end_time;
    start_time =  omp_get_wtime();
    int check_result = find_k_biggest_serial(arr, size, k);
    end_time =  omp_get_wtime();
    present_result(0, size, k, check_result, end_time - start_time);
    return check_result;
}


void copy_array(int * arr_source, int * arr_dest, int size){
    for(int i = 0; i < size; i++){
        arr_dest[i] = arr_source[i];
    }
}
int main() {
    srand(time(NULL));
    int n = 24; // size of the array
    int k = 3; // find the kth biggest number
    int threads = 8;
    int * main_arr;
    int * working_arr;
    int single_result, parallel_result, check_result;


    main_arr = allocate_array(n);
    fill_random(main_arr, n);
    omp_set_num_threads(threads);

    working_arr = allocate_array(n);
    copy_array(main_arr, working_arr, n);
    single_result = calc_single(working_arr, n, k);
    free(working_arr);
    working_arr = allocate_array(n);
    copy_array(main_arr, working_arr, n);
    parallel_result = calc_parallel(working_arr, n, k);
    free(working_arr);

    working_arr = allocate_array(n);
    copy_array(main_arr, working_arr, n);
    check_result = calc_check(working_arr, n, k);
    free(working_arr);


    free(main_arr);
    if (single_result != check_result || parallel_result != check_result)
    {
        return -1;
    }
    

    return 0;
}
