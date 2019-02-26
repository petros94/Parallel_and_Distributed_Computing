#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
struct timeval startwtime, endwtime;
double seq_time;

int number_1=0;
int number_2=0;
int N;          // data array size
int *a;         // data array to be sorted
int TH;

int MAX_THREAD_A=0;
 int MAX_THREAD_B=0;
const int ASCENDING  = 1;
const int DESCENDING = 0;
omp_lock_t lock1;
omp_lock_t lock2;

struct Array{
    int lo;
    int cnt;
    int dir;
};
void init(void);
void print(void);
void sort(void);
void test(void);
void exchange(int i, int j);
void compare(int i, int j, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void recBitonicSortfunc(int lo, int cnt, int dir);
void omprecBitonicSort(int lo, int cnt, int dir);
void ompbitonicMerge(int lo, int cnt, int dir);

/** the main program **/ 
int main(int argc, char **argv) {
    omp_init_lock(&lock1);
    omp_init_lock(&lock2);

    if (argc != 3) {
        printf("Usage: %s arg1 arg2\n  where n=2^arg1 is problem size (power of two) and th=2^arg2 the number of threads\n",
               argv[0]);
        exit(1);
    }
    TH = 1<<atoi(argv[2]);
    N = 1<<atoi(argv[1]);
    if (N<16000){
        MAX_THREAD_A=1;
        MAX_THREAD_B=MAX_THREAD_A;
    }
    else{
        MAX_THREAD_A=TH;
        MAX_THREAD_B=MAX_THREAD_A;
    }
  a = (int *) malloc(N * sizeof(int));

  init();
    //print();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  //print();
   test();
  printf("%f\n", seq_time );
   // printf("number1: %d\n", number_1);
   // printf("number2: %d\n", number_2);
    omp_destroy_lock(&lock1);
    omp_destroy_lock(&lock2);
}

/** -------------- SUB-PROCEDURES  ----------------- **/ 

/** procedure test() : verify sort results **/
void test() {
  int pass = 1;
  int i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1] <= a[i]);
  }

  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}


/** procedure init() : initialize array "a" with data **/
void init() {
  int i;
  for (i = 0; i < N; i++) {
    a[i] = rand() % N; // (N - i);
  }
}

/** procedure  print() : print array elements **/
void print() {
  int i;
  for (i = 0; i < N; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}


/** INLINE procedure exchange() : pair swap **/
 void exchange(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}



/** procedure compare() 
   The parameter dir indicates the sorting direction, ASCENDING 
   or DESCENDING; if (a[i] > a[j]) agrees with the direction, 
   then a[i] and a[j] are interchanged.
**/
 void compare(int i, int j, int dir) {
  if (dir==(a[i]>a[j])) 
    exchange(i,j);
}




/** Procedure bitonicMerge() 
   It recursively sorts a bitonic sequence in ascending order, 
   if dir = ASCENDING, and in descending order otherwise. 
   The sequence to be sorted starts at index position lo,
   the parameter cbt is the number of elements to be sorted. 
 **/

 void bitonicMerge(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    int i;
    for (i=lo; i<lo+k; i++)
      compare(i, i+k, dir);
    bitonicMerge(lo, k, dir);
    bitonicMerge(lo+k, k, dir);
  }
}



/** function recBitonicSort()
    first produces a bitonic sequence by recursively sorting 
    its two halves in opposite sorting orders, and then
    calls bitonicMerge to make them in the same order 
 **/
void recBitonicSortfunc(int lo, int cnt, int dir) {
    int k=cnt/2;
  if (cnt>1) {
      recBitonicSortfunc(lo, k, ASCENDING);
      recBitonicSortfunc(lo+k, k, DESCENDING);
      bitonicMerge(lo, cnt, dir);
  }
}
void ompbitonicMerge(int lo, int cnt, int dir){
    if (cnt>1) {
        int k=cnt/2;
        if (number_2<MAX_THREAD_B-number_1){
            omp_set_lock(&lock2);
            number_2=number_2+2;
            omp_unset_lock(&lock2);
            int i;
            for (i=lo; i<lo+k; i++)
                compare(i, i+k, dir);
            #pragma omp parallel
            {
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        ompbitonicMerge(lo, k, dir);
                    }
                    
                    
                    #pragma omp section
                    {
                        ompbitonicMerge(lo+k, k, dir);
                    }
                    
                }
                #pragma omp barrier
            }
            omp_set_lock(&lock2);
            number_2=number_2-2;
            omp_unset_lock(&lock2);
      }
        else {
            int i;
            for (i=lo; i<lo+k; i++)
                compare(i, i+k, dir);
            bitonicMerge(lo, k, dir);
            bitonicMerge(lo+k, k, dir);
        }
    }
}

void omprecBitonicSort(int lo, int cnt, int dir){
    int k=cnt/2;
    if (cnt>1){
        if (number_1<MAX_THREAD_A){
            omp_set_lock(&lock1);
            number_1=number_1+2;
            omp_unset_lock(&lock1);
            #pragma omp parallel
            {
                int i;
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        omprecBitonicSort(lo, k, ASCENDING);
                    }
                    
                    
                    #pragma omp section
                    {
                        omprecBitonicSort(lo+k, k, DESCENDING);
                    }
                    
                }
                #pragma omp barrier
                #pragma omp master
                {
                    for (i=lo; i<lo+k; i++)
                        compare(i, i+k, dir);
                }
                #pragma omp barrier
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        ompbitonicMerge(lo, k, dir);
                    }
                    
                    
                    #pragma omp section
                    {
                        ompbitonicMerge(lo+k, k, dir);
                    }
                    
                }
                #pragma omp barrier
            }
            omp_set_lock(&lock1);
            number_1=number_1-2;
            omp_unset_lock(&lock1);
        }
        else recBitonicSortfunc(lo, cnt, dir);
    }
}



/** function sort() 
   Caller of recBitonicSort for sorting the entire array of length N 
   in ASCENDING order
 **/
void sort() {
    if (MAX_THREAD_A==1) recBitonicSortfunc(0, N, ASCENDING);
    else {
        omprecBitonicSort(0,N,ASCENDING);
    }
}




