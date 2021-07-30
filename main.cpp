/*
*  This file is part of Christian's OpenMP software lab 
*
*  Copyright (C) 2016 by Christian Terboven <terboven@itc.rwth-aachen.de>
*  Copyright (C) 2016 by Jonas Hahnfeld <hahnfeld@itc.rwth-aachen.de>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>

#include <omp.h>
#include <numa.h>

#include <iostream>
#include <algorithm>

#include <cstdlib>
#include <cstdio>

#include <cmath>
#include <ctime>
#include <cstring>

#define LEFT true
#define RIGHT false

/**
  * tuple type used in heap
  * first: value of element by which it is sorted
  * second: index of origin sub-array where value in first was taken from
  */
struct pair
{
	long first;
	long second;
}pair;

/**
  * heap struct
  * heap: array of heap values
  * capacity: maximum number of elements in heap
  * size: number of elements currently in heap
  */
struct minHeap
{
	struct pair *heap;
	int capacity;
	int size;
};


/**
  * helper routine: check if array is sorted correctly
  */
bool isSorted(int ref[], int data[], const size_t size){
	std::sort(ref, ref + size);
	for (size_t idx = 0; idx < size; ++idx){
		if (ref[idx] != data[idx]) {
			return false;
		}
	}
	return true;
}


/**
 * Binary search, but returning index of searched element instead of Boolean
 * array: array to be searched
 * begin: left delimiter of interval in array to be searched (inclusive)
 * end: right delimiter of interval in array to be searched (exclusive)
 * l_or_r: indicates whether to search for left or right interval end
 * For elements contained in the interval, returns the leftmost occurence if l_or_r is LEFT,
 * 		successor of rightmost occurence if l_or_r is RIGHT
 * For elements not contained in the interval, returns the leftmost element contained if l_or_r is LEFT,
 * 		successor of rightmost element contained if l_or_r is RIGHT
 */
long bin_find(int *array, long begin, long end, int search_item, bool l_or_r){
	long lb = begin-1;
	long rb = end;
	long mid = 0;
	while(rb - lb > 2){
		mid = (lb + rb - 1)/2;
		if(array[mid] == search_item){
			if(l_or_r == LEFT){
				while(array[mid-1] == search_item && mid > begin){
					mid--;
				}
			}
			else{
				while(array[mid+1] == search_item && mid < end){
                    mid++;
                }

			}
			if(l_or_r){
				return mid;
			}
			else{
				return mid+1;
			}
		}
		else if(array[mid] < search_item){
			lb = mid;
		}
		else{
			rb = mid + 1;
		}
	}
	return rb-1;
}

//Hilfsmethoden für den Heap
int getParent(int index) { return (index - 1) / 2; }
int getLeftChild(int index) { return 2 * index + 1; }
int getRightChild(int index) { return 2 * index + 2; }
void swap(struct pair *element1, struct pair *element2){
	struct pair element = *element1;
	*element1 = *element2;
	*element2 = element;
}

void heapify(struct minHeap *in, int index){
	int leftChild = getLeftChild(index);
	int rightChild = getRightChild(index);
	int rightIndex = index;
	if (leftChild < in->size && in->heap[leftChild].first < in->heap[index].first){
		rightIndex = leftChild;
	}
	if (rightChild < in->size && in->heap[rightChild].first < in->heap[rightIndex].first){
		rightIndex = rightChild;
	}
	if (rightIndex != index){
		swap(&in->heap[index], &in->heap[rightIndex]);
		heapify(in, rightIndex);
	}
}

struct minHeap* create_heap(int cap){
	struct minHeap *out = (struct minHeap*) malloc(sizeof(struct minHeap));
	out->heap =(struct pair*) malloc(cap*sizeof(pair));
	out->capacity = cap;
	out->size = 0;
	return out;
}

void delete_heap(struct minHeap* in){
	free(in->heap);
	free(in);
}

bool isEmpty_heap(struct minHeap *in){
	if (in->size > 0){
		return false;
	}
	return true;
}

struct pair extractMin(struct minHeap *in){
	if (in->size == 1){
		in->size--;
		return in->heap[0];
	}
	// Speichert den minimalen Wert und entfernt ihm von Heap
	struct pair root = in->heap[0];
	in->heap[0] = in->heap[in->size - 1];
	in->size--;
	heapify(in, 0);
	return root;
}

bool insert(struct minHeap *in, long neu1, long neu2){
	if (in->size == in->capacity){
		return false;
	}
	in->size++;
	int index = in->size - 1;
	in->heap[index].first = neu1;
	in->heap[index].second = neu2;

	while (index != 0 && in->heap[getParent(index)].first > in->heap[index].first){
		swap(&in->heap[index], &in->heap[getParent(index)]);
		index = getParent(index);
	}

	return true;
}


/**
  * parallel merge step using heap
  * Merges all sorted sub-arrays as specified by sub_starts and sub_ends and writes them to the output array
  * sub_starts: containts first index of corresponding interval in each sorted sub-array
  * sub_starts: contains last index of corresponding interval in each sorted sub-array
  * outBegin: index of out from which to start writing results
  */
void MsMergeParallel(int *out, int *in, long *sub_starts, long *sub_ends, long outBegin, int num_threads) {
	struct minHeap *sub_mins = create_heap(num_threads+1);

	for (int i = 0; i < num_threads; i++){
		if(sub_starts[i] < sub_ends[i]){
			insert(sub_mins, in[sub_starts[i]], i);
			sub_starts[i]++;
		}
	}

	struct pair p;
	int n;

	while(! isEmpty_heap(sub_mins)){
		p = extractMin(sub_mins);

		out[outBegin] = (int) p.first;
		outBegin++;

		n = p.second;

		if(sub_starts[n] < sub_ends[n]){
			insert(sub_mins, in[sub_starts[n]], n);
			sub_starts[n]++;
		}
	}

	delete_heap(sub_mins);
}

/**
  * sequential merge step (straight-forward implementation)
  */
void MsMergeSequential(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
	long left = begin1;
	long right = begin2;

	long idx = outBegin;

	while (left < end1 && right < end2) {
		if (in[left] <= in[right]) {
			out[idx] = in[left];
			left++;
		} else {
			out[idx] = in[right];
			right++;
		}
		idx++;
	}

	while (left < end1) {
		out[idx] = in[left];
		left++, idx++;
	}

	while (right < end2) {
		out[idx] = in[right];
		right++, idx++;
	}
}


/**
  * sequential MergeSort
  */
void MsSequential(int *array, int *tmp, bool inplace, long begin, long end) {
	if (begin < (end - 1)) {
		const long half = (begin + end) / 2;
		MsSequential(array, tmp, !inplace, begin, half);
		MsSequential(array, tmp, !inplace, half, end);
		if (inplace) {
			MsMergeSequential(array, tmp, begin, half, half, end, begin);
		} else {
			MsMergeSequential(tmp, array, begin, half, half, end, begin);
		}
	} else if (!inplace) {
		tmp[begin] = array[begin];
	}
}

/**
  * parallel MergeSort
  * Around 1st parallel loop: split array into num_threads sequential parts, sorting each one with the sequential mergesort implementation
  * Around 2nd and 3rd parallel loop: determine minimum and maximum number contained in array
  * Around 4th parallel loop: split entire interval of numbers into num_threads sub-intervals of equal size
  * 	and save starts and ends of sub-intervals in each sub-array in sub_starts and sub_ends
  * 	and total number of items in each sub-interval in sum_items
  * Around 5th parallel loop: compute where each thread should start writing to output array,
  * 	then merge sub-arrays' parts in corresponding sub_interval via MsMergeParallel
  * 6th parallel loop: copy results from tmp to output array
  */
void MsParallel(int *array, int *tmp, long size, int num_threads){
	long sub_size = ceil((double)size / (double)num_threads);
	#pragma omp parallel for num_threads(num_threads)
	for(int i = 0; i < num_threads; i++){
		if(i != num_threads - 1){
			MsSequential(array, tmp, true, i*sub_size, (i+1)*sub_size);
		}
		else{
			MsSequential(array, tmp, true, i*sub_size, size);
		}
	}

	int min = array[0];
	#pragma omp parallel for reduction(min: min) num_threads(num_threads)
	for(long i = 0; i < num_threads; i++){
		if(array[i*sub_size] < min){
			min = array[i*sub_size];
		}
	}

	int max = array[size-1];
	#pragma omp parallel for reduction(max: max) num_threads(num_threads)
	for(long i = 0; i < num_threads-1; i++){
		if(array[(i+1)*sub_size -1] > max){
			max = array[(i+1)*sub_size -1];
		}
	}

	int interval_size = max - min + 1;
	int sub_interval = ceil((double)interval_size / (double)num_threads);

	long sum_items[num_threads];
	long sub_starts[num_threads][num_threads];
	long sub_ends[num_threads][num_threads];
	#pragma omp parallel for num_threads(num_threads)
	for(int i = 0; i < num_threads; i++){
		sum_items[i] = 0;
		long interval_start = min + i*sub_interval;
		long interval_end = min + (i+1)*sub_interval - 1;	// inclusive interval end
		if(i == num_threads - 1){
			interval_end = max;
		}
		for(int j = 0; j < num_threads-1; j++){
			sub_starts[i][j] = bin_find(array, j*sub_size, (j+1)*sub_size+1, interval_start, LEFT);
			sub_ends[i][j] = bin_find(array, j*sub_size, (j+1)*sub_size+1, interval_end, RIGHT);
			sum_items[i] += sub_ends[i][j] - sub_starts[i][j];
		}
		int j = num_threads-1;
        sub_starts[i][j] = bin_find(array, j*sub_size, size+1, interval_start, LEFT);
        sub_ends[i][j] = bin_find(array, j*sub_size, size+1, interval_end, RIGHT);
        sum_items[i] += sub_ends[i][j] - sub_starts[i][j];
	}

	#pragma omp parallel for num_threads(num_threads)
	for(int i = 0; i < num_threads; i++){
		long write_start = 0;
		for(int j = 0; j < i; j++){
			write_start += sum_items[j];
		}
		MsMergeParallel(tmp, array, sub_starts[i], sub_ends[i], write_start, num_threads);
	}

	#pragma omp parallel for simd num_threads(num_threads)
	for(int i = 0; i < size; i++){
		array[i] = tmp[i];
	}
}


/**
  * MergeSort
  * Chooses between sequential and parallel mergesort implementation depending on number of threads
  * MsParallel might not work with more than sqrt(size) threads,
  * 	number of threads is reduced to this number if necessary and a warning is given
  */
void Mergesort(int *array, int *tmp, const size_t size) {
	int num_threads = omp_get_max_threads();
	if(num_threads*num_threads > size){
		num_threads = (int)sqrt((double)size);
		printf("\nWARNING: sorting %d elements with %d threads leads to a bad block size of less than %d. %d threads will be used instead.\n", size, omp_get_max_threads(), omp_get_max_threads(), num_threads);
	}

	if(num_threads > 1){
		MsParallel(array, tmp, size, num_threads);
	}
	else{
		MsSequential(array, tmp, true, 0, size);
	}
}


/**
  * @brief program entry point
  */
int main(int argc, char* argv[]) {
	// variables to measure the elapsed time
	struct timeval t1, t2;
	double etime;

	// expect one command line arguments: array size
	if (argc != 2) {
		printf("Usage: MergeSort.exe <array size> \n");
		printf("\n");
		return EXIT_FAILURE;
	}
	else {
		const size_t stSize = strtol(argv[1], NULL, 10);
		int *data = (int*) malloc(stSize * sizeof(int));
		int *tmp = (int*) malloc(stSize * sizeof(int));
		int *ref = (int*) malloc(stSize * sizeof(int));

		printf("Initialization...\n");

		srand(95);
		#pragma omp parallel for
		for (size_t idx = 0; idx < stSize; ++idx){
			data[idx] = (int) (stSize * (double(rand()) / RAND_MAX));
		}
		std::copy(data, data + stSize, ref);

		double dSize = (stSize * sizeof(int)) / 1024 / 1024;
		printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);

		gettimeofday(&t1, NULL);
		Mergesort(data, tmp, stSize);
		gettimeofday(&t2, NULL);
		etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
		etime = etime / 1000;

		printf("done, took %f sec. Verification...", etime);
		if (isSorted(ref, data, stSize)) {
			printf(" successful.\n");
		}
		else {
			printf(" FAILED.\n");
		}

		free(data);
		free(tmp);
		free(ref);
	}

	return EXIT_SUCCESS;
}
