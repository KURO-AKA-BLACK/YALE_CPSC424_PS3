#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "timing.h"
#include <unistd.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int main(int argc, char **argv)
{
  double wall_clock_start, cpu_start, wall_clock_end, cpu_end, total_compute_time = 0, total_comm_time = 0, compute_start, compute_end, comm_start, comm_end;
  double ld_clock_start, ld_clock_end;
  int N, i, run, rank, size,b=0;
  double *A, *B, *C, *Ctrue;
  long sizeAB, sizeC;
  int sizes[4] = {1000, 2000, 4000, 8000};
  char files[4][70] = {"/home/cpsc424_ahs3/shared/assignments/assignment3/part2/C-1000.dat",
                       "/home/cpsc424_ahs3/shared/assignments/assignment3/part2/C-2000.dat",
                       "/home/cpsc424_ahs3/shared/assignments/assignment3/part2/C-4000.dat",
                       "/home/cpsc424_ahs3/shared/assignments/assignment3/part2/C-8000.dat"};

  double wctime, Fnorm;

  FILE *fptr;

  
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (run = 0; run < 4; run++){    
    N = sizes[run];
    int lines_per_block = N/size;
    double* temp_C = malloc(N*N * sizeof(double));
    double* send_buf = malloc(((N/size) * (N-N/size+1+N)/2) * sizeof(double));
    double* recv_buf = malloc(((N/size) * (N-N/size+1+N)/2) * sizeof(double));
    double* perma_A = malloc(((N/size) * (N-N/size+1+N)/2) * sizeof(double));
    double* perma_C = malloc((lines_per_block*lines_per_block*size) * sizeof(double));
    int blocks_for_thread[size];
    int displs[size];
    int so_far = 0;
    for (int j = 0; j < size; j++){
      blocks_for_thread[j] = (lines_per_block * (j+1)) * ((lines_per_block * (j+1))+1)/2 - so_far;
      so_far = so_far + blocks_for_thread[j];
      if (j == 0){
        displs[j] = 0;
      }
      else{
        displs[j] = displs[j-1]+blocks_for_thread[j-1];
      }
    }
    
    if (rank == 0){
      N = sizes[run];
      // to make such this header only get to print once
      if (b == 0){
        printf("Matrix multiplication times:\n    N      TIME (secs)    F-norm of Error\n  -----   -------------  -----------------\n");
        b = 1;
      }
      sizeAB = N * (N + 1) / 2; //Only enough space for the nonzero portions of the matrices
      sizeC = N * N;            // All of C will be nonzero, in general!
      
      A = (double *)calloc(sizeAB, sizeof(double));
      B = (double *)calloc(sizeAB, sizeof(double));
      C = (double *)calloc(sizeC, sizeof(double));
      
      srand(12345); // Use a standard seed value for reproducibility
      
      // This assumes A is stored by rows, and B is stored by columns. Other storage schemes are permitted
      for (i = 0; i < sizeAB; i++)
        A[i] = ((double)rand() / (double)RAND_MAX);
      for (i = 0; i < sizeAB; i++)
        B[i] = ((double)rand() / (double)RAND_MAX);
    }
    timing(&wall_clock_start, &cpu_start);  
    comm_start = MPI_Wtime();
    MPI_Scatterv(A, blocks_for_thread, displs, MPI_DOUBLE, perma_A, blocks_for_thread[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, blocks_for_thread, displs, MPI_DOUBLE, recv_buf, blocks_for_thread[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);  
    comm_end = MPI_Wtime();
    total_comm_time += comm_end - comm_start;
    
    for (int loop_num = 0; loop_num < size; loop_num ++){
      if (rank == 0){
        int recv_size = 0;
        compute_start = MPI_Wtime();
        // the if-else are used because recv_size and elements within MIN has different expression between when loop_num is 0 and when loop_num is not 0
        if (loop_num == 0){
          timing(&ld_clock_start, &cpu_start);
          int iC = 0;
          recv_size = blocks_for_thread[0];
          memcpy(send_buf, recv_buf, recv_size * sizeof(double));
          for (int i=0; i<lines_per_block; i++) {
            int iA = i*(i+1)/2; // initializes row pointer in A
            for (int j=0; j<lines_per_block; j++,iC++) {
              int iB = j*(j+1)/2; 
              perma_C[iC] = 0.;
              for (int k=0; k<=MIN(i,j); k++){
                perma_C[iC] += perma_A[iA+k] * recv_buf[iB+k]; 
              }
            }
          }
        }
        else{
          recv_size = (lines_per_block * (size - loop_num+1)) * ((lines_per_block * (size - loop_num+1))+1)/2 - (lines_per_block * (size-1-loop_num+1)) * ((lines_per_block * (size-1-loop_num+1))+1)/2;
          memcpy(send_buf, recv_buf, recv_size * sizeof(double));
          int iC = (size-loop_num)*lines_per_block*lines_per_block;
          for (int i=0; i<lines_per_block; i++) {
            int iA = i*(i+1)/2;
            int cnt = 0;
            for (int j=0; j<lines_per_block; j++,iC++) {
              int iB;
              if (j > 0){
                iB = j * ((size-loop_num)*lines_per_block+1) + cnt;
                cnt = cnt + j;
              }
              else{
                iB = 0;
              }
              perma_C[iC] = 0.;
              for (int k=0; k<=MIN(i,lines_per_block*(size-loop_num)+1+j); k++){
                perma_C[iC] += perma_A[iA+k] * recv_buf[iB+k]; 
              }
            }
          }
        }
        // Calculate the size of incoming data
        int temp = (lines_per_block * (size - loop_num)) * ((lines_per_block * (size - loop_num))+1)/2 - (lines_per_block * (size-1-loop_num)) * ((lines_per_block * (size-1-loop_num))+1)/2;
        compute_end = MPI_Wtime();
        total_compute_time += compute_end - compute_start;
        comm_start = MPI_Wtime();
        // No need to send and receive data after the computation in the last run complete.
        if (loop_num < size - 1){
          MPI_Recv(recv_buf, temp, MPI_DOUBLE, size - 1, 2, MPI_COMM_WORLD, &status);
          MPI_Send(send_buf, recv_size, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
        }
        comm_end = MPI_Wtime();
        total_comm_time += comm_end - comm_start;
        if (loop_num == size - 1){
          timing(&ld_clock_end, &cpu_end);
        }
      }
      
      else{
        if (loop_num == 0){
          timing(&ld_clock_start, &cpu_start);
        }
        compute_start = MPI_Wtime();
        int recv_size = 0;
        // determine the size of incoming column block
        if (rank-loop_num >= 0){
          recv_size = ((lines_per_block * ((rank-loop_num)+ 1)) * ((lines_per_block * ((rank-loop_num) + 1))+1)/2) - ((lines_per_block * (rank-loop_num)) * ((lines_per_block * ((rank-loop_num)))+1)/2);
        }
        else{
          recv_size = ((lines_per_block * (size+(rank-loop_num)+ 1)) * ((lines_per_block * (size+(rank-loop_num) + 1))+1)/2) - ((lines_per_block * (size+rank-loop_num)) * ((lines_per_block * ((size+rank-loop_num)))+1)/2);
        }
        // do job with perma_A and recv_buf
        int t = rank - loop_num;
        if (t < 0){
          t = size+t;
        }
        int iC = t*lines_per_block*lines_per_block;
        int count = 0;
        for (int i=0; i<lines_per_block; i++) {
          int iA;
          if (i == 0){
            iA = 0;
          }
          else{
            iA = i * (rank*lines_per_block + 1) + count;
            count = count + i;
          }
          int cnt = 0;
          for (int j=0; j<lines_per_block; j++,iC++) {
            int iB;
            if (j == 0){
              iB = 0;
            }
            else{
              iB = j * (t*lines_per_block+1) + cnt;
              cnt = cnt + j ;
            }
            perma_C[iC] = 0.;
            for (int k=0; k<=MIN(lines_per_block*rank+i,lines_per_block*t+j); k++){
              perma_C[iC] += perma_A[iA+k] * recv_buf[iB+k]; 
            }
          }
        }
        compute_end = MPI_Wtime();
        total_compute_time += compute_end - compute_start;
        comm_start = MPI_Wtime();
        if (loop_num < size - 1){
          int target = rank + 1;
          if (target == size) target = 0;
          MPI_Send(recv_buf, recv_size, MPI_DOUBLE, target, 2, MPI_COMM_WORLD);
          // calculate the size of incoming data based on the current loop_num
          if (rank-loop_num-1 >= 0){
            recv_size = ((lines_per_block * ((rank-loop_num-1)+ 1)) * ((lines_per_block * ((rank-loop_num-1) + 1))+1)/2) - ((lines_per_block * (rank-loop_num-1)) * ((lines_per_block * ((rank-loop_num-1)))+1)/2);
          }
          else{
            recv_size = ((lines_per_block * (size+(rank-loop_num-1)+ 1)) * ((lines_per_block * (size+(rank-loop_num-1) + 1))+1)/2) - ((lines_per_block * (size+rank-loop_num-1)) * ((lines_per_block * ((size+rank-loop_num-1)))+1)/2);
          }
          MPI_Recv(recv_buf, recv_size, MPI_DOUBLE, rank - 1, 2, MPI_COMM_WORLD, &status);
        }
        comm_end = MPI_Wtime();
        total_comm_time += comm_end - comm_start;
        if (loop_num == size - 1){
          timing(&ld_clock_end, &cpu_end);
        }
      }
    }
    comm_start = MPI_Wtime();
    MPI_Gather(perma_C, lines_per_block*lines_per_block*size, MPI_DOUBLE, temp_C, lines_per_block*lines_per_block*size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    comm_end = MPI_Wtime();
    total_comm_time += comm_end - comm_start;
    timing(&wall_clock_end, &cpu_end);
    // Remainder of the code checks the result against the correct answer (read into Ctrue)
    // reorder the received data
    if (rank == 0){          
      int count = 0;
      for (int f = 0; f < size; f++){
        for (int x = 0; x < lines_per_block; x++){
          for (int i = 0; i < size; i++){
            for (int j = 0; j < lines_per_block; j++){
              C[count] = temp_C[ lines_per_block*N*f +x*lines_per_block + lines_per_block*lines_per_block*i + j];
              count ++;
            }
          }
        }
      }
      
      Ctrue = (double *)calloc(sizeC, sizeof(double));
      
      fptr = fopen(files[run], "rb");
      fread(Ctrue, sizeof(double), sizeC, fptr);
      fclose(fptr);
      
      //Compute the Frobenius norm of Ctrue-C
      Fnorm = 0.;
      for (i = 0; i < sizeC; i++)
        
        Fnorm += (Ctrue[i] - C[i]) * (Ctrue[i] - C[i]);
      Fnorm = sqrt(Fnorm);
      
      // Print a table row
      printf("  %5d    %9.4f  %18.12f\n", N, wall_clock_end-wall_clock_start, Fnorm);
      
      free(Ctrue);

      free(A);
      free(B);
      free(C);
    }
//    printf("rank %d work time: %f, computation time: %f, communication time: %f\n ", rank, ld_clock_end-ld_clock_start, total_compute_time, total_comm_time);
//    MPI_Barrier(MPI_COMM_WORLD);
    free(temp_C);
    free(perma_C);
    free(perma_A);
    free(send_buf);
    free(recv_buf);
  }
  MPI_Finalize();
}
