#include "ransac_2Dline.h"
#include <sys/time.h>

int main(int argc, char **argv) {

	char *path = NULL;
	int N = 0;

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

	float data[2*N];

	FILE * pdata;
	pdata = fopen(path, "r");

	if(pdata == NULL)
	{
		perror("can't open file");
		exit(1);
	}

	char line[10000];
	int i = 0;
	while(!feof(pdata))
	{
		if(fgets(line, 10000, pdata))
		{
			char *tok = strtok(line, ",");
			while(tok!=NULL)
			{
				data[i++] = atof(tok);
				tok = strtok(NULL, ",");
			}
		}
	}

//	for(int j=0; j<2*N; j++)
//		printf("%.2f  ",data[j]);

	float bestModel[3];
	int bestInliers;
	struct timeval t1, t2;
	double elapsedTime;

	// start timer
    gettimeofday(&t1, NULL);

	ransac_2Dline(data, N, N, 149, 0.3, bestModel, &bestInliers, 1);

	// stop timer
    gettimeofday(&t2, NULL);

	// compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

	if(bestInliers>0)
	printf("\n<<<<<<<<<< RANSAC END >>>>>>>>>>>>>\n model %.2f*x + %.2f*y + %.2f = 0\n time= %f ms\n inliers = %d\n\n", bestModel[0], 
bestModel[1], bestModel[2], elapsedTime, bestInliers);

	return 0;
}
