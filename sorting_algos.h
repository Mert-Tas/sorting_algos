#ifndef SORTING_ALGOS_H
#define SORTING_ALGOS_H

#ifdef BUILD_DLL
    #define DLL_EXPORT __declspec(dllexport)
#else
    #define DLL_EXPORT __declspec(dllimport)
#endif

#ifdef __cplusplus
extern "C"
{
#endif

class Sorting_Algos
{
    public:
        Sorting_Algos();
        virtual ~Sorting_Algos();

        int DLL_EXPORT bubbleSort(int *a, int n);
        int DLL_EXPORT selectionSort(int *a, int n);
        int DLL_EXPORT insertionSort(int *a, int n);
        void DLL_EXPORT insertionSortRec(int *a, int n);
        int DLL_EXPORT shellSort(int *a, int n);
        void DLL_EXPORT quickSort(int *a, int left, int rigth);
        void DLL_EXPORT mergeSort(int *a, int left, int right);
        void DLL_EXPORT heapSort(int *a, int n);
        void DLL_EXPORT bucketSort(float *a, int n);
        void DLL_EXPORT radixSort(int *a, int n);

    private:
        void XORswap(int &a, int &b);
        void swapElementsByRef(int &a, int &b);
        // Partitioning function for quick sort
        int partitionArray(int *a, int left, int rigth);
        int medianOfThree(int *a, int left, int right);
        // Merge function for merge sort a given array and indices left, middle, right
        void mergeArrays(int *a, int l, int m, int r);
        // Make an array max heap
        void maxHeapify(int *a, int n, int i);
        // Count sort. Helper sort for radix sort
        void countSort(int *a, int n, int us);
        // Gets the maximum element in an array. Helper for radix sort
        int getMaxArray(int *a, int n);

};  // End of class Sorting_Algos

#ifdef __cplusplus
}
#endif

#endif // SORTING_ALGOS_H
