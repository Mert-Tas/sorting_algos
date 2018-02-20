#include "sorting_algos.h"
#include <algorithm>   // max
#include <limits>      // limits in C++

using namespace std;

DLL_EXPORT Sorting_Algos::Sorting_Algos()
{
    //ctor
}

DLL_EXPORT Sorting_Algos::~Sorting_Algos()
{
    //dtor
}

void Sorting_Algos::swapElementsByRef(int &a, int &b)
{
    int temp = a;
    a = b;
    b = temp;
}

void Sorting_Algos::XORswap(int &a, int &b)
{
    a = a ^ b; // Store both a and b in hybrid a. XORing elements stores them
    b = b ^ a; // Cancel b value from hybrid a, leaving a. [(b^b = 0) ^ a] = a
    a = a ^ b; // Cancel a value from hybrid a, leaving b. [(a^a = 0) ^ b] = b
}

// Starting on the left, compare adjacent items and keep "bubbling"
// the larger one to the right (it’s in its final place).
// Bubble sort the remaining N -1 items.
// Best Case : O(n)   Worst Case : O(n^2)
int DLL_EXPORT Sorting_Algos::bubbleSort(int *a, int n)
{
    int cnt = 0;

    for (int i = n; i >= 1; --i)
    {
        // Check if the inner loop is already sorted so that sorting is finished
        bool swapped = false;
        for (int j = 1; j < i; ++j)
        {
            if (a[j-1] > a[j])
            {
                XORswap(a[j-1], a[j]);
                swapped = true;
            }
        }

        cnt ++;
        if (!swapped)
            break;
    }

    return cnt;
}

// Scan all items and find the smallest. (comparing starting from first)
// Swap it into position as the first item.
// Repeat the selection sort on the remaining N-1 items.
// Best Case / Worst Case : O(n^2)
int DLL_EXPORT Sorting_Algos::selectionSort(int* a, int n)
{
    int cnt = 0;

    for (int i = 0; i < n; ++i)
    {
        int minInd = i;
        bool minFound = false;

        for (int j = i+1; j < n; ++j)
        {
            if (a[j] < a[minInd])
            {
                minInd = j;
                minFound = true;
            }
        }

        if (minFound)
        {
            XORswap(a[i], a[minInd]);
            cnt ++;
        }
    }

    return cnt;
}

// Start with a sorted list of 1 element on the left,
// and N-1 unsorted items on the right. Take the first
// unsorted item (element #2) and insert it into the sorted list,
// moving elements as necessary.
// We now have a sorted list of size 2, and N -2 unsorted elements.
// Repeat for all elements.
// Best Case : O(n)   Worst Case : O(n^2)
int DLL_EXPORT Sorting_Algos::insertionSort(int *a, int n)
{
    int key, j;
    int cnt = 0;
    // Assume first element is sorted so start from the second element.
    for (int i = 1; i < n; ++i)
    {
        key = a[i];
        j = i-1;

        // Inner loop: sorts elements in the left of index i  (element 0 and 1)
        while (j >= 0 && a[j] > key)
        {
            a[j+1] = a[j];
            j --;
        }
        // Move to the new element a[i]
        a[j+1] = key;
        cnt++;
    }

    return cnt;
}

// Recursive insertion sort
void DLL_EXPORT Sorting_Algos::insertionSortRec(int* a, int n)
{
    // Base case: if array size is one or smaller, return
    if (n < 1)
        return;

    // Sort first n-1 items
    insertionSortRec(a, n-1);

    // Insert last element at its correct position
    int last = a[n-1];
    int j = n-2;

    // Move elements of arr[0..i-1], that are
    // greater than key, to one position ahead
    // of their current position
    while (j >= 0 && a[j] > last)
    {
        a[j+1] = a[j];
        j --;
    }
    a[j+1] = last;
}

// ShellSort is mainly a variation of Insertion Sort.
// In insertion sort, we move elements only one position ahead.
// When an element has to be moved far ahead, many movements are involved.
// The idea of shellSort is to allow exchange of far items.
// In shellSort, we make the array h-sorted for a large value of h.
// We keep reducing the value of h until it becomes 1.
// An array is said to be h-sorted if all sublists of every h'th element is sorted.
// Complexity : O(n^2)
int DLL_EXPORT Sorting_Algos::shellSort(int *a, int n)
{
    int cnt = 0;

    // Start with a big gap, then reduce the gap
    for (int gap = n/2; gap > 0; gap /= 2)
    {
        // Do a gapped insertion sort for this gap size.
        // The first gap elements a[0..gap-1] are already in gapped order
        // keep adding one more element until the entire array is
        // gap sorted
        for (int i = gap; i < n; ++i)
        {
            // add a[i] to the elements that have been gap sorted
            // save a[i] in temp and make a hole at position i
            int temp = a[i];

            // shift earlier gap-sorted elements up until the correct
            // location for a[i] is found
            int j;
            for (j = i; j >= gap && a[j - gap] > temp; j -= gap)
                a[j] = a[j - gap];

            // put temp (the original a[i]) in its correct location
            a[j] = temp;
            cnt ++;
        }
    }
    return cnt;
}

// The main function that implements QuickSort
// arr[]: Array to be sorted,
// left:  Starting index,
// right: Ending index
// Best Case : O(n logn)   Worst Case : O(n^2)
// Together with heap sort, quick sort is not stable.
void DLL_EXPORT Sorting_Algos::quickSort(int *a, int left, int right)
{
    if (left < right)
    {
        // p is the partitioning index a[p] is at the right place
        int p = partitionArray(a, left, right);

        // Recursively sort the left and right sub arrays
        quickSort(a, left, p-1);
        quickSort(a, p+1, right);
    }
}

int Sorting_Algos::medianOfThree(int *a, int left, int right)
{
    int middle = (left + (right - left)) / 2;
    int median  = max(a[middle], max(a[left], a[right]));

    // Pick the middle of three
    if (median == a[middle])
        return max(a[left], a[right]);
    if (median == a[left])
        return max(a[middle], a[right]);
    if (median == a[right])
        return max(a[left], a[middle]);

    return a[middle];
}

// This function takes last element as pivot, places
// the pivot element at its correct position in sorted
// array, and places all smaller (smaller than pivot)
// to left of pivot and all greater elements to right of pivot
int Sorting_Algos::partitionArray(int *a, int left, int right)
{
    // Element to be positioned at the right position
    int pivot = a[right];

    // Index of smaller element
    int i = left - 1;

    for (int j = left; j <= right - 1; ++j)
    {
        // Compare current element with pivot
        if (a[j] <= pivot)
        {
            // Increment index of smaller element
            ++i;
            swapElementsByRef(a[i], a[j]);
        }
    }
    // The pivot is often swapped to the front,
    // so it is out of the way during the pivoting.
    // Afterwards, it is swapped into place
    // (with a pivot item that is less than or equal to it,
    // so the pivot is preserved).
    swapElementsByRef(a[i+1], a[right]);

    return i+1;
}

// Merges two sub arrays of a[].
// First sub array is a[l..m]
// Second sub array is a[m+1..r]
void Sorting_Algos::mergeArrays(int *a, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temp arrays
    int L[n1], R[n2];

    // Copy data to temp arrays L and R
    for (int i = 0; i < n1; ++i)
        L[i] = a[l + i];
    for (int j = 0; j < n2; ++j)
        R[j] = a[m + 1 + j];

    // Merge the temp arrays back into a[1..r]
    i = 0;  // Initial index of first sub array
    j = 0;  // Initial index of second sub array
    k = l;  // Initial index of merged sub array is left
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            a[k] = L[i];
            ++i;
        }
        else
        {
            a[k] = R[j];
            ++j;
        }
        ++k;
    }

    // Copy the remaining elements of L[], if there are any
    while (i < n1)
        a[k++] = L[i++];

    // Copy the remaining elements of R[], if there are any
    while (j < n2)
        a[k++] = R[j++];
}

// It divides input array in two halves, calls itself for the two halves and
// then merges the two sorted halves.
// Useful for sorting linked lists.
// Time Complexity: O(n logn) for all cases. log2n for recursion and n times for sorting
// Space Complexity: O(n). It is not a sort in place algorithm.
void DLL_EXPORT Sorting_Algos::mergeSort(int *a, int left, int right)
{
    if (left < right)
    {
        // Find the middle element. (r - l) to avoid overflow for large l and r
        int middle = left + (right - left) / 2;

        // Sort first and second halves
        mergeSort(a, left, middle);
        mergeSort(a, middle + 1, right);

        mergeArrays(a, left, middle, right);
    }
}

void Sorting_Algos::maxHeapify(int *a, int n, int i)
{
    // Initialize the nodes
    // Root of the heap (binary tree) should be the largest element
    int root  = i;
    int left  = (i << 1) + 1;     //  left = 2*i + 1
    int right = (i + 1) << 1;     // right = 2*i + 2 = 2*(i+1)

    // Compare the root to make it the largest element.
    // If left is bigger than the root
    if (left < n && a[left] > a[root])
        root = left;

    // If right is bigger than the root
    if (right < n && a[right] > a[root])
        root = right;

    // After comparisons if root is not the largest element
    // swap root
    if (root != i)
    {
        swapElementsByRef(a[root], a[i]);
        // Recursively heapify the affected sub tree
        maxHeapify(a, n, root);
    }
}

// Similar to quick sort, its typical implementation is not stable, but can be made stable
// In-place sorting algorithm
// Time complexity: O(n logn)
// Create and build heap is O(logn)
// Quick sort and merge sort are better in practice.
void DLL_EXPORT Sorting_Algos::heapSort(int *a, int n)
{
    // Build heap (rearrange array)
    // Start from bottommost and rightmost internal mode and heapify all
    // internal modes in bottom up way
    for (int i = n/2 - 1; i >= 0; --i)
        maxHeapify(a, n, i);

    // One by one extract an element from heap
    for (int i = n-1; i >= 0; --i)
    {
        // Move current root (largest element) to end of array
        swapElementsByRef(a[0], a[i]);

        // Call max heapify on the reduced heap
        maxHeapify(a, i, 0);
    }
}

// Array elements are between 0.0f and 1.0f
void DLL_EXPORT Sorting_Algos::bucketSort(float *a, int n)
{
    // 1) Create n empty buckets O(n)
    vector <float> buc[n];

    // 2) Put array elements in different buckets O(n)
    for (int i = 0; i < n; ++i)
    {
        // Index in bucket
        int bucind = n * a[i];
        buc[bucind].push_back(a[i]);
    }

    // 3) Sort individual buckets  O(n)
    for (int i = 0; i < n; ++i)
        sort(buc[i].begin(), buc[i].end());

    // 4) Concatenate all buckets into a[]   O(n) ? as there will be n items in all buckets.
    int index = 0;
    for (int i = 0; i < n; ++i)
        for (size_t j = 0; j < buc[i].size(); ++j)
            a[index++] = buc[i][j];
}

// Helper sort function for radix sort
void Sorting_Algos::countSort(int *a, int n, int exp)
{
    int output[n];
    // There are 10 digits
    int digit[10] = {0};

    // Stores each digit frequency in the corresponding index of digit array
    // Ex: [3, 7, 23, 111]
    // exp = 1 (ones digit)  : digit[0] = 2 {3, 7}, digit[1] = 1 {111}, digit[3] = 1 {23}
    // exp = 10 (tens digit) : digit[1] = 1 {111} , digit[2] = 1 {23}
    for (int i = 0; i < n; ++i)
    {
        ++ digit[(a[i] / exp) % 10];
    }

    // Modify the count array such that each element at each index
    // stores the sum of previous counts.
    // The modified count array indicates the position of each object
    // in the output sequence. Loop to 10
    for (int i = 1; i < 10; ++i)
    {
        digit[i] += digit[i-1];
    }

    // Output each object from the input sequence followed by
    // decreasing its count by 1. Reverse loop, normal loop does not works.
    for (int i = n - 1; i >= 0; --i)
    {
        output[digit[(a[i] / exp) % 10] - 1] = a[i];
        digit[(a[i] / exp) % 10] --;
    }

    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (int i = 0; i < n; ++i)
    {
        a[i] = output[i];
    }
}

int Sorting_Algos::getMaxArray(int *a, int n)
{
    int maxEl = a[0];

    for (int i = 1; i < n; ++i)
    {
        if (a[i] > maxEl)
            maxEl = a[i];
    }
    return maxEl;
}

// Radix Sort sorts numbers digit by digit starting from least
// significant digit to most significant digit.
// Radix sort uses counting sort as a subroutine to sort
// Complexity : Let there be d digits in input integers.
// Radix Sort takes O(d*(n+b)) time where b is the base for
// representing numbers, for example, b is 10 for decimal system.
// If we have log2n bits for every digit, the running time of Radix
// appears to be better than Quick Sort for a wide range of input numbers.
void DLL_EXPORT Sorting_Algos::radixSort(int *a, int n)
{
    // Find the maximum number to know the number of digits
    int maxEl = getMaxArray(a, n);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (int exp = 1; maxEl/exp > 0; exp *= 10)
    {
        countSort(a, n, exp);
    }
}
