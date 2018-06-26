#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "heat.h"
#include "helper.h"
//#include <stdbool.h>

double calculate_total_heatx(double *h_new, int x_per_th,int y_per_th)
{
    double heat = 0.0; // total heat in system
    for (int j = 1; j < y_per_th + 1; ++j)
    {
      for (int i = 1; i < x_per_th + 1; ++i)
      {
        heat += h_new[map(i, j, x_per_th+2)];
      }
    }
    return heat;
}

void printarr_x(double *a, int x_per_th, int  y_per_th, int rank) {
  // does nothing right now, should record each "frame" as image
  char name[20];
  sprintf(name, "heat_%d.svg", rank);
  FILE *fp = fopen(name, "w");
  const int size = 5;

  fprintf(fp, "<html>\n<body>\n<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">");

  fprintf(fp, "\n<rect x=\"0\" y=\"0\" width=\"%i\" height=\"%i\" style=\"stroke-width:1;fill:rgb(0,0,0);stroke:rgb(0,0,0)\"/>", size*x_per_th, size*y_per_th);
  for(int j=1; j<y_per_th+1; ++j)
    for(int i=1; i<x_per_th+1; ++i) {
      int rgb = (a[map(i,j,x_per_th+2)] > 0) ? rgb = (int)round(255.0*a[map(i,j,x_per_th+2)]) : 0.0;
      if(rgb>255) rgb=255;
      if(rgb) fprintf(fp, "\n<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"stroke-width:1;fill:rgb(%i,0,0);stroke:rgb(%i,0,0)\"/>", size*(i-1), size*(j-1), size, size, rgb, rgb);
    }
  fprintf(fp, "</svg>\n</body>\n</html>");

  fclose(fp);
}

/*
* calculates a 4 size list of bool, which indicates, if the process with "rank" is in the boundaries of stencil.
* mapping of the list is (top,right,bottom,left). For example, if the process is the most right of the stencil
* then it calculates (. , true, . , .).
*/
void boundary_binary(int* boundaries, int px, int py, int rank){
   int reminder = rank % px;
   int divide = rank/px;
   if(rank < px){
     boundaries[0] =  1; // top
   }
   else{
     boundaries[0] = 0;
   }
   if(reminder == (px-1)){
     boundaries[1] =  1; // right
   }
   else{
     boundaries[1] = 0;
   }
   if(divide == (py-1)){
     boundaries[2] =  1; // bottom
   }
   else{
     boundaries[2] = 0;
   }
   if(reminder == 0){
     boundaries[3] =  1; // left
   }
   else{
     boundaries[3] = 0;
   }

}

void traverse_inner(int x_per_th, int y_per_th, double* h_new, double* h_old ){
  int n = x_per_th +2;
  for(int j=2; j<y_per_th; j++){ // row count
    for(int i=2; i<x_per_th; i++){ // column count //after traversing each inner point in a row, go to next row. so the cache acc. is utilized.
      h_new[map(i, j ,n)] = h_old[map(i, j, n)] / 2.0 + (h_old[map(i - 1, j, n)] + h_old[map(i + 1, j, n)] + h_old[map(i, j - 1, n)] + h_old[map(i, j + 1, n)]) / 8.0 ;
    }
  }
}

void traverse_col(int x_per_th, int y_per_th, double* h_new, double* h_old,
                  int row_start_ind, int row_count, int col_ind ){
  int row_end = row_start_ind + row_count;
  int n = x_per_th +2;
  for(int j=row_start_ind; j<row_end; j++){ // row count
    h_new[map(col_ind, j ,n)] = h_old[map(col_ind, j, n)] / 2.0 + (h_old[map(col_ind - 1, j, n)] + h_old[map(col_ind + 1, j, n)] + h_old[map(col_ind, j - 1, n)] + h_old[map(col_ind, j + 1, n)]) / 8.0;
  }
}
void traverse_col_nonbound(int x_per_th, int y_per_th, double* h_new, double* h_old,
                  int row_start_ind, int row_count, int col_ind, MPI_Request* req ){ //traverse non boundary columns(not boundary in whole matrix. so this column might be first column of a right matrix.)
  int row_end = row_start_ind + row_count;
  int n = x_per_th +2;
  for(int j=row_start_ind; j<row_end; j++){ // row count
    MPI_Wait(&req[j-row_start_ind],MPI_STATUS_IGNORE);
    h_new[map(col_ind, j ,n)] = h_old[map(col_ind, j, n)] / 2.0 + (h_old[map(col_ind - 1, j, n)] + h_old[map(col_ind + 1, j, n)] + h_old[map(col_ind, j - 1, n)] + h_old[map(col_ind, j + 1, n)]) / 8.0;
  }
}
void traverse_row(int x_per_th, int y_per_th, double* h_new, double* h_old,
                  int col_start_ind, int col_count, int row_ind ){
  int col_end = col_start_ind + col_count;
  int n = x_per_th +2;
  for(int i=col_start_ind; i<col_end; i++){ // row count
    h_new[map(i, row_ind ,n)] = h_old[map(i, row_ind, n)] / 2.0 + (h_old[map(i - 1, row_ind, n)] + h_old[map(i + 1, row_ind, n)] + h_old[map(i, row_ind - 1, n)] + h_old[map(i, row_ind + 1, n)]) / 8.0;
  }
}

double jacobi(double *h_new, double *h_old, int niters, int energy_intensity,
              int n, int iter_energy,  const int nsources, int sources[nsources][2],
               int rank, int size, int px, int py, MPI_Comm comm, int output){

    int x_per_th = n / px;
    int y_per_th = n / py;

    int the_row_count = x_per_th + 2;
    h_old = (double *)calloc(1, (the_row_count) * (y_per_th + 2) * sizeof(double)); // extended with halos of width 1
    h_new = (double *)calloc(1, (the_row_count) * (y_per_th + 2) * sizeof(double)); // extended with halos of width 1

    //double* temp_col1 = (double *)calloc(1, (y_per_th) * sizeof(double));
    //double* temp_col2 = (double *)calloc(1, (y_per_th) * sizeof(double));

    int* boundaries= (int*) malloc(sizeof(int) * 4);
    boundary_binary(boundaries, px, py, rank);

    // to send top, need 1 request, since horizontal array
    // to send right, need y_per_th requests, each is a number and not horizontal
    // to send bottom, need 1 request, horizontal array
    // to send left, need y_per_th requests, each is a number and not horizontal


    double *tmp;
    int top_req_ind = 0;
    int right_req_ind = top_req_ind + (!boundaries[0])*(1);
    int bottom_req_ind = right_req_ind + (!boundaries[1])*(y_per_th);
    int left_req_ind = bottom_req_ind + (!boundaries[2])*(1);


    const int num_requests_per_iter = (!boundaries[0])*(1)+ (!boundaries[1])*(y_per_th) + (!boundaries[2])*(1) + (!boundaries[3])*(y_per_th);
    for(int iter=0; iter<niters; ++iter){

      MPI_Request req[num_requests_per_iter*2];
      //printf("rank:%d,iter:%d\n", rank,iter);
      // MPI_Barrier(comm);
      // if(rank==0){
      //   printf("iter:%d\n",iter );
      // }
      if(boundaries[1]!=1){ // if the thread is not in the right boundary then send the last column to the right...
        int base_send = iter*y_per_th;
        for(int i=1; i<y_per_th+1; i++){
          //printf("!!!rank:%d right_req_ind:%d num_requests_per_iter:%d\n",rank ,right_req_ind,num_requests_per_iter);
          MPI_Isend(&h_old[map(x_per_th, i ,the_row_count)], 1,MPI_DOUBLE, rank+1,base_send+i,comm,&req[i-1+right_req_ind]);
          MPI_Irecv(&h_old[map(x_per_th+1, i ,the_row_count)], 1,MPI_DOUBLE, rank+1,base_send+i,comm,&req[i-1+right_req_ind+num_requests_per_iter]);

        }
      }

      if(boundaries[2] != 1){ // if the thread is not in the bottom boundary then send the last row to the bottom...
        MPI_Isend(&h_old[map(1, y_per_th ,the_row_count)], x_per_th,MPI_DOUBLE, rank+px,iter,comm,&req[bottom_req_ind]);
        MPI_Irecv(&h_old[map(1, y_per_th+1 ,the_row_count)], x_per_th ,MPI_DOUBLE, rank+px,iter,comm,&req[bottom_req_ind+num_requests_per_iter]);
      }
      if(boundaries[0] != 1){ // if the thread is not in the top boundary then send the first row to the upper one...
        MPI_Isend(&h_old[map(1, 1 ,the_row_count)], x_per_th,MPI_DOUBLE, rank-px,iter,comm,&req[top_req_ind]);
        MPI_Irecv(&h_old[map(1, 0 ,the_row_count)], x_per_th ,MPI_DOUBLE, rank-px,iter,comm,&req[top_req_ind+num_requests_per_iter]);
      }
      if(boundaries[3]!=1){ // if the thread is not in the left boundary then send the first column to the left...
        int base_send = iter*y_per_th;
        for(int i=1; i<y_per_th+1; i++){
          MPI_Isend(&h_old[map(1, i ,the_row_count)], 1,MPI_DOUBLE, rank-1,base_send+i,comm,&req[i-1+left_req_ind]);
          MPI_Irecv(&h_old[map(0, i ,the_row_count)], 1,MPI_DOUBLE, rank-1,base_send+i,comm,&req[i-1+left_req_ind+num_requests_per_iter]);
        }
      }

      traverse_inner(x_per_th, y_per_th, h_new, h_old ); // first traverse the part where no dependency needed.

      if(boundaries[0]!=1){
        MPI_Wait(&req[top_req_ind+num_requests_per_iter],MPI_STATUS_IGNORE);
      }
      traverse_row(x_per_th, y_per_th, h_new, h_old, 2, x_per_th-2, 1 ); // traverse first row, without left and right boundaries.
      if(boundaries[2]!=1){
        MPI_Wait(&req[bottom_req_ind+num_requests_per_iter],MPI_STATUS_IGNORE);
      }
      traverse_row(x_per_th, y_per_th, h_new, h_old, 2, x_per_th-2, y_per_th ); // traverse last row, without left and right boundaries.


      if(boundaries[3]!=1){ // if not most left, then it is a nonboundary
        traverse_col_nonbound(x_per_th,y_per_th, h_new, h_old,1, y_per_th, 1, &req[left_req_ind+num_requests_per_iter] );// traverse first(left) column
      }
      else{// if most left, then it is a boundary column
        traverse_col(x_per_th,y_per_th, h_new, h_old,1, y_per_th, 1 );// traverse first(left) column      }
      }
      if(boundaries[1]!=1){ // if not most right, then it is a nonboundary
        traverse_col_nonbound(x_per_th,y_per_th, h_new, h_old,1, y_per_th, x_per_th, &req[right_req_ind+num_requests_per_iter] );// traverse first(left) column
      }
      else{
        traverse_col(x_per_th,y_per_th, h_new, h_old,1, y_per_th, x_per_th );// traverse first(left) column
      }

      if (iter < iter_energy){
        for (int i = 0; i < nsources; ++i){
          int rank_i = (sources[i][0]-1) / x_per_th;
          int rank_j = (sources[i][1]-1) / y_per_th;
          if(rank_j*px + rank_i == rank){
            if(iter==0){
              printf("rank:%d source:%d, 0:%d, 1:%d\n",rank,i,sources[i][0],sources[i][1] );
              printf("follow up: i:%d, j:%d\n",((sources[i][0]-1)%x_per_th) + 1,((sources[i][1]-1)%y_per_th) + 1 );
            }
            h_new[map(((sources[i][0]-1)%x_per_th) + 1, ((sources[i][1]-1)%y_per_th) + 1, x_per_th+2)] += energy_intensity; // heat rate
          }
        }
      }
      tmp = h_new; // swap arrays
      h_new = h_old;
      h_old = tmp;

      MPI_Waitall(num_requests_per_iter, &req[0], MPI_STATUS_IGNORE);
    }

    if (output) printarr_x(h_old, x_per_th,y_per_th, rank);
    return calculate_total_heatx(h_old, x_per_th,y_per_th);
}
