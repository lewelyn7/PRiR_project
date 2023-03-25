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

int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = low - 1;

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] > pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}
struct bounds{
    int down;
    int up;
};
int quickselect(int arr[], int size, int low, int high, int k, int threads) {
    int * p_indexes = allocate_array(threads);
    struct bounds * bounds_arr = malloc(threads * sizeof(struct bounds));
    for(int i=0; i<threads; i++){
		bounds_arr[i].down = low + i * size / threads;
		bounds_arr[i].up = i < threads - 1 ? low + (i+1) * size / threads - 1 : high;
    }
    for(int i=0; i<threads; i++){
        p_indexes[i] = partition(arr, low, high);
    }
    if (low == high) {
        return arr[low];
    }
    int pi = partition(arr, low, high);
    int rank = pi - low + 1;

    if (k == rank) {
        return arr[pi];
    } else if (k < rank) {
        return quickselect(arr, low, pi - 1, k);
    } else {
        return quickselect(arr, pi + 1, high, k - rank);
    }
}

int* allocate_array(int n) {
    int* arr = (int*) malloc(n * sizeof(int));
    return arr;
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
    #pragma omp parallel shared(parallel_result)
    {
        #pragma omp single nowait
        {
            parallel_result = quickselect(arr, 0, size-1, k);
        }
    }
    end_time =  omp_get_wtime();
    present_result(2, size, k, parallel_result, end_time - start_time);
    return parallel_result;
}

int calc_single(int * arr, int size, int k){
    double start_time, end_time;
    start_time =  omp_get_wtime();
    int simple_result = quickselect(arr, 0, size-1, k);
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
    int n = 10000000; // size of the array
    int k = 50; // find the kth biggest number
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
