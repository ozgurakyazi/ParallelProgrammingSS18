#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <pthread.h>

#include "mandelbrot_set.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))


struct man_args
{
	int x_resolution;
	int y_start;
	int y_end;
	int max_iter;
	double view_x0;
	double view_x1;
	double view_y0;
	double view_y1;
	double x_stepsize;
	double y_stepsize;
	int palette_shift;
	void* img;//[x_resolution][3];
};

/*
void update_image(int x_resolution,unsigned char (*img)[x_resolution][3],char color[3] ){

}
*/
void* traverse(void* args){

	struct man_args* arg= (struct man_args*) args;
	double y;
	double x;

	unsigned char (*img)[arg->x_resolution][3] = (unsigned char (*)[arg->x_resolution][3])(arg->img);

	complex double Z;
	complex double C;

	int k;

	for (int i = arg->y_start; i < arg->y_end; i++){
		for (int j = 0; j < arg->x_resolution; j++){
			y = arg->view_y1 - i * arg->y_stepsize;
			x = arg->view_x0 + j * arg->x_stepsize;

			Z = 0 + 0 * I;
			C = x + y * I;

			k = 0;

			do{
				Z = Z * Z + C;
				k++;
			} while (cabs(Z) < 2 && k < arg->max_iter);

			if (k == arg->max_iter){
				memcpy(img[i][j], "\0\0\0", 3);
				//printf("%s haaa \n", img[i][j]);
			}
			else{
				int index = (k + arg->palette_shift)
				            % (sizeof(colors) / sizeof(colors[0]));
				memcpy(img[i][j], colors[index], 3);
				//printf("%s i:%d , j:%d \n", img[i][j],i,j);
			}
		}
	}
	return NULL;
}
void mandelbrot_draw(int x_resolution, int y_resolution, int max_iter,
	                double view_x0, double view_x1, double view_y0, double view_y1,
	                double x_stepsize, double y_stepsize,
	                int palette_shift, unsigned char (*image)[x_resolution][3],
						 int num_threads) {

		pthread_t *threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));

		struct man_args* args = (struct man_args*) malloc(num_threads * sizeof(struct man_args) );

		int section_size = y_resolution / num_threads;
		//int resid = (y_resolution % num_threads != 0) ? 1 : 0;

		for(int i=0; i < num_threads ; i++){
				args[i].x_resolution = x_resolution;
				args[i].y_start = i*section_size;
				args[i].y_end = num_threads!=i+1 ? ((i+1)*section_size)  : y_resolution;
				// printf("###############################\n" );
				// printf("%d\n",args[i].y_start );
				// printf("%d\n",args[i].y_end );
				// printf("###############################\n" );
				args[i].max_iter = max_iter;
				args[i].view_x0 = view_x0;
				args[i].view_x1 = view_x1;
				args[i].view_y0 = view_y0;
				args[i].view_y1 = view_y1;
				args[i].x_stepsize = x_stepsize;
				args[i].y_stepsize = y_stepsize;
				args[i].palette_shift = palette_shift;
				args[i].img = image;
				pthread_create(&threads[i], NULL, traverse, args+i);
		}

		for(int i=0; i<num_threads; i++){
			pthread_join(threads[i],NULL);
		}

		free(threads);
		free(args);
}
