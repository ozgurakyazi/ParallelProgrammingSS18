#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "helper.h"

void reverse(char *str, int strlen)
{
        // parallelize this function and make sure to call reverse_str()
        // on each processor to reverse the substring.

        int np, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //MPI_Status s[np-1];
    char *my_str = str;
    if(rank!=0){
      my_str = (char *) malloc(sizeof(char) * (strlen+1));
    }

    int new_strlen = strlen;
    int resid = strlen % np;
    int n_padding = (np - resid)%np;
    if(resid != 0){
      new_strlen += n_padding;
    }
    int chunk_size = new_strlen/np;
    int base_ind = chunk_size-n_padding;
    if(n_padding > chunk_size){
      chunk_size = (strlen-resid)/np;
      base_ind = resid+chunk_size;
    }
    MPI_Bcast(my_str, strlen, MPI_CHAR, 0, MPI_COMM_WORLD);

    int new_rank = np - rank -1;
    int start_ind = new_rank * chunk_size; // start index is inclusive
    //char *result = NULL;
    //char *to_return = NULL;
    if (rank == 0){

        //printf("the chunk_size is:%d\n", chunk_size);
        //printf("the n_padding is:%d\n", n_padding);
        //result = (char *) malloc(sizeof(char) * (strlen));
        //reverse_str(&my_str[start_ind], base_ind);
        if(np!=1){
          for(int i=0;i<base_ind; i++){
            str[i] = str[start_ind+base_ind-i-1];
          }
        }
        else{
          reverse_str(&my_str[start_ind], base_ind);
          memcpy(str,&my_str[start_ind],(base_ind)*(sizeof(char)) );
        }
        //to_return = &my_str[start_ind-n_padding];
        // MPI_Gather(&str[start_ind-n_padding],chunk_size, MPI_CHAR, result
        //  ,chunk_size ,MPI_CHAR, 0, MPI_COMM_WORLD);
        //memcpy(str,&my_str[start_ind],(base_ind)*(sizeof(char)) );
        //free(result);
        int current_index = base_ind;
        for(int i=1; i<np; i++){
          //printf("%d is waiting for %d\n", rank,i);
          MPI_Recv(&str[current_index], chunk_size,MPI_CHAR, i,15,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          //printf("%d is done with %d\n", rank,i);
          current_index+= chunk_size;

        }
        //memcpy(&str[base_ind],&result[base_ind],(chunk_size*(np-1))*(sizeof(char)) );
        //printf("the result is:%d\n", (strlen-resid));
    }
    else{
      reverse_str(&my_str[start_ind], chunk_size);
      //to_return = &my_str[start_ind];
      //printf("%d is waiting and my result is %s \n", rank,&my_str[start_ind]);
      MPI_Send(&my_str[start_ind], chunk_size,MPI_CHAR, 0,15,MPI_COMM_WORLD);
      //printf("%d is done\n", rank);
      //MPI_Gather(&my_str[start_ind],chunk_size, MPI_CHAR, result
      //  ,chunk_size ,MPI_CHAR, 0, MPI_COMM_WORLD);
      //free(my_str);
    }


    //MPI_Gather(to_return,chunk_size, MPI_CHAR, result
    // ,chunk_size ,MPI_CHAR, 0, MPI_COMM_WORLD);

    // if(rank==0){
    //   memcpy(str,&result[n_padding],strlen*(sizeof(char)) );
    //   free(result);
    // }else{
    //
    // }

}
