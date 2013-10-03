#include "ransac_2Dline.h"
#include "svd.h"
#include <sys/time.h>

#define BSIZE 10000

int main(int argc, char **argv) {

	char *path = NULL;
	int N = 0;
	int i;

	if(argc!=3)
	{
		perror("USE: ransac path_to_dataset n\n");
		exit(1);
	}
	else
	{
		path = argv[1];
		N = atoi(argv[2]);
	}

	printf("path=%s ; N=%d\n",path,N);

	double **data = malloc(N * sizeof(double *));
	if(data == NULL) { perror("out of memory\n"); exit(0); }
	
	for(i = 0; i < N; i++)
	{
		data[i] = malloc(2 * sizeof(double));
		if(data[i] == NULL) { perror("out of memory\n"); exit(0); }
	}

	FILE * pdata;
	pdata = fopen(path, "r");

	if(pdata == NULL) { perror("can't open file"); exit(1); }

	char buffer[BSIZE];
	int k = 0;
	while(!feof(pdata))
	{
		if(fgets(buffer, BSIZE, pdata))
		{
			char *tok = strtok(buffer, ",");
			while(tok!=NULL)
			{
				printf("k = %d\n",k);
				if(k < N)
					data[k++][0] = atof(tok);
				else
					data[(k++)-N][1] = atof(tok);
				
				tok = strtok(NULL, ",");
			}
		}
	}

	printf("Data =\n"); printMatrix(data,N,2);

	double bestModel[3];
	int bestInliers;
	
	struct timeval t1, t2;
	double elapsedTime;

	// start timer
    gettimeofday(&t1, NULL);

	ransac_2Dline(data, N, (N/2)-1, 1, bestModel, &bestInliers, 1);

	// stop timer
    gettimeofday(&t2, NULL);

	// compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

	if(bestInliers>0)
	printf("\n<<<<<<<<<< RANSAC END >>>>>>>>>>>>>\n model (%lf)*x + (%lf)*y + (%lf) = 0\n time= %f ms\n inliers = %d\n\n", bestModel[0], 
bestModel[1], bestModel[2], elapsedTime, bestInliers);

	return 0;
}
