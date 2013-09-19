#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


void ransac_2Dline(float *data, int n, int nd, int maxT, float threshold,
					float *bestModel, int *bestInliers, int verbose);

int randomSelect(float *sel, int nsel, float *data, int *ndata, int N);

int fitModel_line(float *point, float *model, float threshold);

void estimateModel_line(float *model, float *dataset, int n, int inliers);
