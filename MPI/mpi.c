#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <sys/time.h>

int volume;
int n,m,k;
int processNumber;
int numberOfPoints;
int *processMap;

struct Point{
  int id;
  int boxId;
  int boxCoordinates[3];
  float pointCoordinates[3];
  int closestNeighbourId;
  int neighbourBoxId;
  float distance;
};

struct Box{
  int boxId;
  struct Point *points;
};

void CreateRandomPoint(struct Point *pointer, int id, int n, int m, int k, struct Box *boxes );
float Distance(struct Point *point1, struct Point *point2);
void InitializeGrid(int volume,int *n, int *m, int *k);
void AssignPointsToBoxes(struct Point *pointset, struct Box *box);
void AssignBoxesToProcesses(struct Box *box, struct Box **process);
void InitializeProcess(struct Box **a);
void ArrangeBoxes(struct Point *pointer);
void ShowBoxId(struct Point point);
void ShowPointCoordinates(struct Point point);
int WhichProcess(int boxId);
int WhichBox(int boxId, struct Box *box);
int Convert(int x, int y, int z);
void CreateMap(struct Point *point, int *map);
void CreateMPIpointtype(MPI_Datatype *point_type);
void SendToProcess(struct Box **P_q);
void ReceiveBoxes(struct Box *boxes, int rank, int* counter);
int CountElements(struct Point *point);
struct Point LocalSearch(struct Point *point, struct Box *box, int counterC);
int CheckBounds(struct Point point);
void SendOutbox(struct Point **outbox);
void ReceiveInbox(struct Point **inbox);
int CheckPoint(struct Point *pointReceived, struct Box *boxes);
int cmp(void *a, void *b);

MPI_Datatype point_type;
struct timeval startwtime, endwtime;
double seq_time1, seq_time2;

int main(int argc, char *argv[]){
    numberOfPoints=1<<atoi(argv[1]);
    volume=1<<atoi(argv[2]);
    MPI_Init(&argc, &argv);
    int  numtasks, rank;
    CreateMPIpointtype(&point_type);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    processNumber=numtasks;
    time_t t;
    srand((unsigned) time(&t));   //random number generation, every second
    //numberOfPoints=2*65536;
    //volume=0.5*0.5*0.0625*4096;
    InitializeGrid(volume, &n, &m, &k);
    if (rank==0){          //Initial work is done by first process
      struct Point *q;
      struct Point *c;
      struct Box *boxessQ;
      struct Box *boxessC;
      struct Box **P_q;
      struct Box **P_c;
      P_q=malloc(processNumber*sizeof(struct Box*));
      P_c=malloc(processNumber*sizeof(struct Box*));
      boxessQ=malloc(volume*sizeof(struct Box));
      boxessC=malloc(volume*sizeof(struct Box));
      for (int i=0; i<volume; i++){
        boxessQ[i].points=malloc(4*(numberOfPoints/volume)*sizeof(struct Point));
        boxessQ[i].boxId=i;
        boxessC[i].points=malloc(4*(numberOfPoints/volume)*sizeof(struct Point));
        boxessC[i].boxId=i;
      }
      q=malloc(numberOfPoints*sizeof(struct Point));
      c=malloc(numberOfPoints*sizeof(struct Point));

      for (int i=0; i<numberOfPoints; i++){
          CreateRandomPoint(&q[i],0, n, m, k, boxessQ);
          CreateRandomPoint(&c[i],1, n, m, k, boxessC);
      }
      gettimeofday (&startwtime, NULL);
      InitializeProcess(P_q);
      InitializeProcess(P_c);
      /*for (int i=0; i<numberOfPoints; i++){
          ArrangeBoxes(&q[i]);
          ArrangeBoxes(&c[i]);
      }*/
      //qsort(q, numberOfPoints, sizeof(struct Point), cmp);
      //qsort(c, numberOfPoints, sizeof(struct Point), cmp);
      int counterQ=0;
      int counterC=0;
      for (int i=0; i<volume; i++){
        AssignPointsToBoxes(q, &boxessQ[i]);
        AssignPointsToBoxes(c, &boxessC[i]);
        int p=CountElements(boxessQ[i].points);
        for (int k=0; k<p; k++){
          if (boxessQ[i].points[k].pointCoordinates[0]==0) printf("Error, please run again.\n" );
        }
      }
      AssignBoxesToProcesses(boxessQ, P_q);
      AssignBoxesToProcesses(boxessC, P_c);
      SendToProcess(P_q);
      SendToProcess(P_c);
      MPI_Barrier(MPI_COMM_WORLD);
      gettimeofday (&endwtime, NULL);
      for (int i=1; i<processNumber; i++){
        free(P_q[i]);
        free(P_c[i]);
      }



      /* After initializing random points, placing them to boxes,
      assinging those boxes to processes and sending them over, each process
      now starts finding closest neighbours for all points assigned to it. */
      seq_time1 = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);
      printf("%f\n", seq_time1);
      gettimeofday (&startwtime, NULL);

      int *countBoxesQ;
      int *countBoxesC;

      /* Creating first outbox. Contains points whose neighbours qualify to
      be in other boxes*/

      struct Point **outboxOne;
      outboxOne=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        outboxOne[i]=malloc(128*volume*sizeof(struct Point));
      }
      countBoxesQ= malloc(volume/processNumber*sizeof(int));
      countBoxesC= malloc(volume/processNumber*sizeof(int));

      for (int i=0; i<volume/processNumber; i++){
        countBoxesQ[i]= CountElements(P_q[0][i].points);
        countBoxesC[i]= CountElements(P_c[0][i].points);
      }

/* Section below contains code which is used to fill the outbox. For every
single point, the process checks if it's located close to the boundaries
of the box. If it does, it finds which process is in charge of the nearest
box, and places the point in the outbox addressed to it.*/

      int map[27];
      int counter=0;
      int flag=1;
      for (int m=0; m<processNumber; m++){
        for (int i=0; i<volume/processNumber; i++){
          for (int j=0; j<countBoxesQ[i]; j++){
            flag=1;
            if (CheckBounds(P_q[0][i].points[j])==1){
              CreateMap(&P_q[0][i].points[j], map);
              for (int k=0; k<27; k++){
                if (map[k]>=0 && k!=13){
                  if (WhichProcess(map[k])==m && flag==1) {
                    P_q[0][i].points[j].neighbourBoxId=map[k];
                    outboxOne[m][counter]= P_q[0][i].points[j];
                    counter++;
                    //flag=0;
                  }
                }
              }
            }
          }
        }
        //printf("Process %d sent %d points to %d \n",rank, counter, m );
        counter=0;
      }

      SendOutbox(outboxOne);

/* After sending the outbox, the process now receives points from other processes.
For each point received, the process finds the closest neigbour and puts it in
a new outbox. After this is done for every point received, it sends back the
outbox to the other processes */

      for (int i=0; i<volume/processNumber; i++){
        for (int k=0; k<countBoxesQ[i]; k++){
          struct Point possibleCandidate= LocalSearch(&P_q[0][i].points[k], &P_c[0][i], countBoxesC[i]);
          P_q[0][i].points[k].closestNeighbourId=possibleCandidate.id;
          P_q[0][i].points[k].distance= Distance(&possibleCandidate, &P_q[0][i].points[k]);
        }
      }

      struct Point **inboxOne;
      inboxOne=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        inboxOne[i]=malloc(128*volume*sizeof(struct Point));
      }

      ReceiveInbox(inboxOne);
      //printf("%d %d\n", elements, rank);

      struct Point **outboxTwo;
      outboxTwo=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        outboxTwo[i]=malloc(128*volume*sizeof(struct Point));
      }

      for (int i=0; i<processNumber; i++){
        int elements=CountElements(inboxOne[i]);
        //printf("Process %d received %d points from %d \n",rank , elements,i);
        for (int j=0; j<elements; j++){
          int a=WhichBox(inboxOne[i][j].neighbourBoxId, P_q[0]);
          //printf("%d\n",a);
          struct Point pointToSend=LocalSearch(&inboxOne[i][j], &P_c[0][a], countBoxesC[a]);
          pointToSend.closestNeighbourId=inboxOne[i][j].id;
          pointToSend.boxId=inboxOne[i][j].boxId;
          outboxTwo[i][j]=pointToSend;
          //ShowPointCoordinates(outboxTwo[i][j]);
        }
        //ShowPointCoordinates(outboxTwo[i][elements-1]);
      }

      printf("%d\n",CountElements(inboxOne[0]) );
      SendOutbox(outboxTwo);

      struct Point **inboxTwo;
      inboxTwo=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        inboxTwo[i]=malloc(128*volume*sizeof(struct Point));
      }

      ReceiveInbox(inboxTwo);
      int f=0;
      int d=0;
      printf("%d\n",CountElements(inboxTwo[0]) );
      for (int i=0; i<processNumber; i++){
        int b=CountElements(inboxTwo[i]);
        for (int j=0; j<b; j++){
          //printf("%d\n",inboxTwo[i][j].boxId);
          f=f+CheckPoint(&inboxTwo[i][j], P_q[0]);
        }
        d=d+b;
      }
      printf("Changed %d out of %d points\n", f, d);
      for (int i=0; i<volume/processNumber; i++){
        for (int j=0; j<countBoxesQ[i]; j++){
          if (P_q[0][i].points[j].closestNeighbourId==0) printf("jeronimo\n");
        }
      }

      printf("\n");
      MPI_Barrier(MPI_COMM_WORLD);
      gettimeofday (&endwtime, NULL);
      seq_time2 = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);
      printf("%f\n", seq_time2);
      printf("%f\n", seq_time1);
      MPI_Barrier(MPI_COMM_WORLD);

    }
    else {
      struct Box *boxesQ;
      struct Box *boxesC;
      boxesQ=malloc(volume/processNumber*sizeof(struct Box));
      boxesC=malloc(volume/processNumber*sizeof(struct Box));
      int *countBoxesQ;
      int *countBoxesC;
      int *countOutbox;
      countBoxesQ= malloc(volume/processNumber*sizeof(int));
      countBoxesC= malloc(volume/processNumber*sizeof(int));
      ReceiveBoxes(boxesQ, rank, countBoxesQ);
      ReceiveBoxes(boxesC, rank, countBoxesC);
      MPI_Barrier(MPI_COMM_WORLD);
      /* work starts here */


      struct Point **outboxOne;
      outboxOne=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        outboxOne[i]=malloc(128*volume*sizeof(struct Point));
      }


      countOutbox=malloc(processNumber*sizeof(int));

      int map[27];
      int counter=0;
      int flag=1;
      for (int m=0; m<processNumber; m++){
        for (int i=0; i<volume/processNumber; i++){
          for (int j=0; j<countBoxesQ[i]; j++){
            if (CheckBounds(boxesQ[i].points[j])==1){
              CreateMap(&boxesQ[i].points[j], map);
              for (int k=0; k<27; k++){
                if (map[k]>=0 && k!=13){
                  if (WhichProcess(map[k])==m && flag==1) {
                    boxesQ[i].points[j].neighbourBoxId=map[k];
                    outboxOne[m][counter]= boxesQ[i].points[j];
                    counter++;
                    //flag=0;
                  }
                }
              }
            }
            flag=1;
          }
        }
        //printf("Process %d sent %d points to %d \n",rank, counter, m );
        counter=0;
      }

      SendOutbox(outboxOne);
      for (int i=0; i<volume/processNumber; i++){
        //printf("%d\n", counterQ);
        for (int k=0; k<countBoxesQ[i]; k++){
          struct Point possibleCandidate= LocalSearch(&boxesQ[i].points[k], &boxesC[i], countBoxesC[i]);
        }
      }

      struct Point **inboxOne;
      inboxOne=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        inboxOne[i]=malloc(128*volume*sizeof(struct Point));
      }

      ReceiveInbox(inboxOne);

      struct Point **outboxTwo;
      outboxTwo=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        outboxTwo[i]=malloc(128*volume*sizeof(struct Point));
      }

      for (int i=0; i<processNumber; i++){
        int elements=CountElements(inboxOne[i]);
        //printf("Process %d received %d points from %d \n",rank , elements,i);
        for (int j=0; j<elements; j++){
          int a=WhichBox(inboxOne[i][j].neighbourBoxId, boxesC);
          struct Point pointToSend=LocalSearch(&inboxOne[i][j], &boxesC[a], countBoxesC[a]);
          pointToSend.closestNeighbourId=inboxOne[i][j].id;
          pointToSend.boxId=inboxOne[i][j].boxId;
          outboxTwo[i][j]=pointToSend;
        }
      }

      SendOutbox(outboxTwo);
      struct Point **inboxTwo;
      inboxTwo=malloc(processNumber*sizeof(struct Point*));
      for (int i=0; i<processNumber; i++){
        inboxTwo[i]=malloc(128*volume*sizeof(struct Point));
      }

      ReceiveInbox(inboxTwo);

      for (int i=0; i<processNumber; i++){
        int b=CountElements(inboxTwo[i]);
        for (int j=0; j<b; j++){
          //printf("%d\n",inboxTwo[i][j].boxId);
          int rc=CheckPoint(&inboxTwo[i][j], boxesQ);
        }
      }

      printf("\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }
}

/*Function for assigning random coordinates
to points and placing them to the grid*/

void CreateRandomPoint(struct Point *pointer, int id, int n, int m, int k, struct Box *boxes) {
  for (int i=0; i<3; i++){
    pointer->pointCoordinates[i]=0.0001*(rand() % 10000)+0.0001;
  }

  for (float j=0; j<n; j++){
    if (pointer->pointCoordinates[0]> j/n){
      pointer->boxCoordinates[0]= j;
    }
    else break;
  }

  for (float j=0; j<m; j++){
    if (pointer->pointCoordinates[1]> j/m){
      pointer->boxCoordinates[1]= j;
    }
    else break;
  }

  for (float j=0; j<k; j++){
    if (pointer->pointCoordinates[2]> j/k){
      pointer->boxCoordinates[2]= j;
    }
    else break;
  }
  pointer->id=id;
  ArrangeBoxes(pointer);

}

/* Distance calculation*/

float Distance(struct Point *point1, struct Point *point2){
  float distance= sqrt((point1->pointCoordinates[0]-point2->pointCoordinates[0])*
  (point1->pointCoordinates[0]-point2->pointCoordinates[0])+(point1->pointCoordinates[1]
  -point2->pointCoordinates[1])*(point1->pointCoordinates[1]-point2->pointCoordinates[1])
  + (point1->pointCoordinates[2]-point2->pointCoordinates[2])*(point1->pointCoordinates[2]
  -point2->pointCoordinates[2]));
  return distance;
}

/*Initialize grid as cubic as possible*/

void InitializeGrid(int volume, int *n, int *m, int *k){
    int a=log2(volume);
    *n= pow(2, (int)(a/3));
    *m= pow(2, (int)((a-(int)(a/3))/2));
    *k= pow(2, a-(int)(a/3)-(int)((a-(int)(a/3))/2));
}

/* Arrange Boxes in linear order based on cartesian coordinates.
The order is Z>Y>X and each position is calculated as
A = N*M*z + N*y + z, where N is max(x) and M is max(y) */

void ArrangeBoxes(struct Point *pointer){
    pointer->boxId= n*m*pointer->boxCoordinates[2]+n*pointer->boxCoordinates[1]
    +pointer->boxCoordinates[0];
}

/* Prints box coordinates and number after arrangement, in console */

void ShowBoxId(struct Point point){
  printf("box x y z: %d %d %d\n",point.boxCoordinates[0], point.boxCoordinates[1],
  point.boxCoordinates[2]);
  printf("box id: %d\n", point.boxId);
}

/* Prints Point Coordinates in console */

void ShowPointCoordinates(struct Point point){
  printf("point x y z: %f %f %f\n", point.pointCoordinates[0], point.pointCoordinates[1],
   point.pointCoordinates[2]);
}

void AssignPointsToBoxes(struct Point *pointset, struct Box *box){
  int counter=0;
  for (int i=0; i<numberOfPoints; i++){
      if (pointset[i].boxId==box->boxId){
        box->points[counter]=pointset[i];
        counter++;
      }
  }
}

void AssignBoxesToProcesses(struct Box *box, struct Box **process){
  int *counter;
  counter=malloc(processNumber*sizeof(int));
  for (int i=0; i<processNumber; i++){
    counter[i]=0;
  }
  for (int i=0; i<volume; i++){
    for (int j=0; j<processNumber; j++){
      if ((box[i].boxId< (j+1)*volume/processNumber) &&
    (box[i].boxId>= j*volume/processNumber)){
      process[j][counter[j]]=box[i];
      counter[j]++;
      break;
      }
    }
  }
}

/* Initializes memory space for each process. Initial size of each
process is chosen to be the same for all processes and equal to number
of points divided by number of processes */

void InitializeProcess(struct Box **a){
  for (int i=0; i<processNumber; i++){
    a[i]=malloc((volume/processNumber)*sizeof(struct Box));
  }
}

int WhichProcess(int boxId){
  for (int i=0; i<processNumber; i++){
    if (boxId<(i+1)*volume/processNumber) return i;
  }
  return -1;
}

int Convert(int x, int y, int z){
  if (z<k && z >=0 &&
    y <m && y >=0 &&
    x <n && x >=0){
      return n*m*z+n*y+x;
    }
    else return -1;
}

void CreateMap(struct Point *point, int *map){
  int counter=0;
  for (int i=-1; i<2; i++){
    for (int j=-1; j<2; j++){
      for (int k=-1; k<2; k++){
        map[counter]=Convert(point->boxCoordinates[0]+k, point->boxCoordinates[1]+j, point->boxCoordinates[2]+i);
        counter++;
      }
    }
  }
}

void SendToProcess(struct Box **P_q){
  for (int i=1; i<processNumber; i++){
    for (int j=0; j<volume/processNumber; j++){
      int count= CountElements(P_q[i][j].points);
      MPI_Send(P_q[i][j].points, count, point_type, i, 0, MPI_COMM_WORLD);
      //printf("I sent box %d to process %d\n", j, i);
      MPI_Send(&P_q[i][j].boxId, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      //printf("and and the box id too!\n");
    }
  }
}

void ReceiveBoxes(struct Box *boxes, int rank, int *counter){
  for (int i=0; i<volume/processNumber; i++){
    MPI_Status status;
    int count, boxId;
    struct Point *incomingBuffer;
    incomingBuffer=malloc(100*numberOfPoints*sizeof(struct Point));
    MPI_Recv(incomingBuffer, 100*numberOfPoints, point_type, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, point_type, &count);
    boxes[i].points=malloc(count*sizeof(struct Point));
    for (int j=0; j<count; j++){
      boxes[i].points[j]=incomingBuffer[j];
    }
    counter[i]=count;
    MPI_Recv(&boxId, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    boxes[i].boxId=boxId;
    free(incomingBuffer);
  }
}

void CreateMPIpointtype(MPI_Datatype *point_type){
  int blocks[7]={1,1,3,3,1,1,1};
  MPI_Datatype types[7]={MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_INT, MPI_INT, MPI_FLOAT};
  MPI_Aint displacements[7];
  MPI_Aint intex, floatex;
  MPI_Type_extent(MPI_FLOAT, &floatex);
  MPI_Type_extent(MPI_INT, &intex);
  displacements[0]= (MPI_Aint)(0);
  displacements[1]= intex;
  displacements[2]= 2*intex;
  displacements[3]= 5*intex;
  displacements[4]= 5*intex+3*floatex;
  displacements[5]= 6*intex+3*floatex;
  displacements[6]= 6*intex+4*floatex;
  MPI_Type_struct(7, blocks, displacements, types, point_type);
  MPI_Type_commit(point_type);
}

int CountElements(struct Point *point){
  int counter=0;
  while (1){
    if ((point[counter].pointCoordinates[0]==(float)0) &&
    (point[counter+1].pointCoordinates[0]==(float)0)) {
        break;
      }
      counter++;
  }
  return counter;
}

struct Point LocalSearch(struct Point *point, struct Box *box, int counterC){
  float min=1;
  int id=0;
  for (int j=0; j<counterC; j++){
    if (Distance(&box->points[j], point)<min){
      min=Distance(&box->points[j], point);
      id=j;
    }
  }
  return box->points[id];
}

int WhichBox(int boxId, struct Box *box){
  for (int i=0; i<volume/processNumber; i++){
    if (boxId==box[i].boxId) return i;
  }
  return -1;
}


int CheckBounds(struct Point point){
  if ((point.pointCoordinates[0]< (float)volume/(numberOfPoints*n)+(float)point.boxCoordinates[0]/n) ||
  (point.pointCoordinates[0] > ((float)point.boxCoordinates[0]+1)/n -(float)volume/(numberOfPoints*n))) {
      return 1;
  }
  else if ((point.pointCoordinates[1]< (float)volume/(numberOfPoints*m)+(float)point.boxCoordinates[1]/m) ||
  (point.pointCoordinates[1] > ((float)point.boxCoordinates[1]+1)/m -(float)volume/(numberOfPoints*m))) {
    return 1;
  }
  else if ((point.pointCoordinates[2]< (float)volume/(numberOfPoints*k)+(float)point.boxCoordinates[2]/k) ||
  (point.pointCoordinates[2] > ((float)point.boxCoordinates[2]+1)/k -(float)volume/(numberOfPoints*k))) {
    return 1;
  }
  else return 0;
}

void SendOutbox(struct Point **outbox){
  MPI_Request req;
  for (int i=0; i<processNumber; i++){
    int e=CountElements(outbox[i]);
    MPI_Isend(outbox[i], e, point_type, i, 0, MPI_COMM_WORLD, &req);
  }
}

void ReceiveInbox(struct Point **inbox){
  MPI_Status status;
  for (int i=0; i<processNumber; i++){
    MPI_Recv(inbox[i], 128*volume, point_type, i, 0, MPI_COMM_WORLD, &status);
  }
}


int CheckPoint(struct Point *pointReceived, struct Box *boxes){
  int counter1, counter2;
  int q=WhichBox(pointReceived->boxId, boxes);
  int elements=CountElements(boxes[q].points);
    for (int j=0; j<elements; j++){
      if (pointReceived->closestNeighbourId==boxes[q].points[j].id){
        if (Distance(pointReceived, &boxes[q].points[j]) < boxes[q].points[j].distance){
          //printf("i changed it. Previous distance: %f\n",boxes[q].points[j].distance);
          boxes[q].points[j].closestNeighbourId=pointReceived->id;
          boxes[q].points[j].distance=Distance(pointReceived, &boxes[q].points[j]);
          //printf("now:%f\n", boxes[q].points[j].distance);
          return 1;
        }
        else return 0;

      }
    }
    return 0;

}

int cmp(void *a, void *b){
    struct Point *pointer1=(struct Point *)a;
    struct Point *pointer2=(struct Point *)b;
    return (pointer1->boxId - pointer2->boxId);
}
