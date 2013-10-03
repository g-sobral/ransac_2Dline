#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int ransac_2Dline(double **data, int n, int maxT, double threshold,
					double *bestModel, int *bestInliers, int verbose);

int randomSelect(double **sel, int nsel, double **data, int *ndata);

int fitModel_line(double *point, double *l, double threshold);

void estimateModel_line(double *l, double **P, int n);
