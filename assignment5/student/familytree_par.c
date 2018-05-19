#include "familytree.h"
#include <omp.h>
int limit = 20;

void traverse2(tree *node){
	if(node != NULL){

		//#pragma omp parallel
		//#pragma omp single
		//{

				#pragma omp task
				{
					node->IQ = compute_IQ(node->data);
					genius[node->id] = node->IQ;
					//printf("thread:%d computes iq as:%d\n", omp_get_thread_num(),node->IQ);
				}
				#pragma omp task
				{
					//printf("%dthread for right-start:\n", omp_get_thread_num());
					traverse2(node->right);
				}

				#pragma omp task
				{
					//printf("%dthread for left-start:\n", omp_get_thread_num());
					traverse2(node->left);
				}
			//}
			//#pragma omp taskwait


	}
}

void traverse(tree *node, int numThreads){

	#pragma omp master
	{
		//omp_set_nested(1);
		//omp_set_max_active_levels(3);

	}

		//#pragma omp sections nowait
		#pragma omp parallel shared(node) num_threads(numThreads)
		{
			#pragma omp single nowait
			traverse2(node);
			#pragma omp taskwait



		}




}
