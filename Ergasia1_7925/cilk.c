#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
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
void cilkrecBitonicSort(int lo, int cnt, int dir);
void cilkbitonicMerge(int lo, int cnt, int dir);

/** the main program **/ 
int main(int argc, char **argv) {

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
  printf("%f\n", seq_time);
    //printf("number1: %d\n", number_1);
    //printf("number2: %d\n", number_2);
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
void cilkbitonicMerge(int lo, int cnt, int dir){
    if (cnt>1) {
        int i;
        int k=cnt/2;
        if (number_2<MAX_THREAD_A-number_1){
            number_2=number_2+1;
            for (i=lo; i<lo+k; i++) compare(i, i+k, dir);
            cilk_spawn cilkbitonicMerge(lo, k ,dir);
            cilkbitonicMerge(lo+k , k ,dir);
            cilk_sync;
            number_2=number_2-1;
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

void cilkrecBitonicSort(int lo, int cnt, int dir){
    int k=cnt/2;
    if (cnt>1){
        int i;
        if (number_1<MAX_THREAD_A){
            number_1=number_1+1;
            cilk_spawn cilkrecBitonicSort(lo, k ,ASCENDING);
            cilkrecBitonicSort(lo+k , k ,DESCENDING);
            cilk_sync;
	    number_1=number_1-1;
            number_2=number_2+1;
            for (i=lo; i<lo+k; i++) compare(i, i+k, dir);
            cilk_spawn cilkbitonicMerge(lo, k, dir);
            cilkbitonicMerge(lo+k, k, dir);
            cilk_sync;
            number_2=number_2-1;
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
        cilkrecBitonicSort(0,N,ASCENDING);
    }
}




