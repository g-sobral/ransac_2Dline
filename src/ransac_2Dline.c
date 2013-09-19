#include "ransac_2Dline.h"

void ransac_2Dline(float *data, int n, int nd, int maxT, float threshold,
					float *bestModel, int *bestInliers, int verbose) {

	if(verbose)
		printf("Start RANSAC, n=%d, maxT=%d, t=%.2f\n\n", n, maxT, threshold);

	*bestInliers = 0;

	int T = 1;
	int ndata = nd;
	int	inliers = 0;
	int Tcount = 0;

	float fracInliers = 0;
//	float conSet[2*n];
	float maybeInliers[4];
	float maybeModel[3];
	float point[2];
	float pNoOutliers = 0;
	float p = 0.99;

	// create a copy of the data array
	float dataCpy[2*nd];
	for(int i=0; i<nd; i++)
	{
		dataCpy[i] = data[i];
		dataCpy[nd+i] = data[n+i];
	}

	srand(time(NULL)); // set rand seed

	while(T > Tcount)
	{
		if(verbose)
			printf("\n#%d ITERATION >>>>>>>>>>>>>>>>>>>>>>>\n", Tcount);

		// Select 2 points at random to form a trial model
		if(randomSelect(maybeInliers, 2, dataCpy, &ndata, nd)==-1)
			break;

		if(verbose)
			printf(" selected points: (%.3f, %.3f) and (%.3f, %.3f)\n", 
					maybeInliers[0], maybeInliers[2], maybeInliers[1], 
					maybeInliers[3]);

		// Fit model to the random selection of data points
		estimateModel_line(maybeModel, maybeInliers, 2, 2);

		if(verbose)
			printf(" model:  %.3f*x + %.3f*y + %.3f = 0\n", maybeModel[0], 
					maybeModel[1], maybeModel[2]);

		inliers = 0;
		fracInliers = 0;

		// Evaluate distances between points and model.
		// Given a threshold, create a consensus set with the points
		// that are inliers.
		for(int i=0; i<nd; i++)
		{
			point[0] = data[i];
			point[1] = data[n+i];

			if(fitModel_line(point, maybeModel, threshold))
			{
//				conSet[inliers] = point[0];
//				conSet[n+inliers] = point[1];
				inliers++;
			}
		}

		if(verbose)
			printf(" inliers = %d\n", inliers);

		if(inliers > *bestInliers)	// Largest set of inliers.
		{
			// Record data for this model
			bestModel[0] = maybeModel[0];
			bestModel[1] = maybeModel[1];
			bestModel[2] = maybeModel[2];
			*bestInliers = inliers;

			// Reestimate T, the number of trials to ensure we pick,
			// with probability p, a data set free of outliers.
			fracInliers = (float)inliers/nd;
			pNoOutliers = 1 - pow(fracInliers, 2);
			T = log(1-p)/log(pNoOutliers);
		}

		Tcount++;
		if(Tcount > maxT)
			break;
	}

	if(bestInliers==0)
	{
		if(verbose)
			printf("\n### ERROR: ransac was unable to find a useful \
solution.\n");
		exit(1);
	}

	if(verbose)
		printf("\n>>>>>>>>>>> RANSAC END\n bestModel %.3f*x + %.3f*y + %.3f \
= 0\n bestInliers=%d\n", bestModel[0], bestModel[1], 
				bestModel[2], *bestInliers);
}

int randomSelect(float *sel, int nsel, float *data, int *ndata, int N) {

    int r = 0;
    int k = *ndata;

	if(nsel > *ndata)
	{
		printf("randomSelect: unable to select %d points from dataset[%d]\n", 
				nsel, *ndata);
		return -1;
	}

	for(int i=0; i<nsel; i++, k--)
    {
        r = rand()%(k);

        sel[i] = data[r];
        sel[nsel+i] = data[N+r];

        data[r] = data[k-1];
        data[N+r] = data[N+k-1];
    }

	*ndata = k;
	printf("ndata = %d\n", *ndata);

	return 0;
}

int fitModel_line(float *point, float *model, float threshold) {
    // Estimate distance between point and model
    // d = abs(a*x + b*y +c)/sqrt(a^2 + b^2)

    float d=0;

    d = fabs(model[0]*point[0] + model[1]*point[1] + model[2])/ \
        sqrt(pow(model[0], 2) + pow(model[1], 2));

    if(d<=threshold)
        return 1;
    else
        return 0;
}

void estimateModel_line(float *model, float *dataset, int n, int inliers) {
    // Estimate the best fitting line by linear regression.
    // equation of line a*x + b*y + c = 0

    float sumX=0, sumY=0, sumXY=0, sumXX=0;

    for(int i=0; i<inliers; i++)
    {
        sumX += dataset[i];
        sumY += dataset[n+i];
        sumXY += dataset[i]*dataset[n+i];
        sumXX += dataset[i]*dataset[i];
    }

    model[0] = (sumXY - (sumX*sumY)/n)/(sumXX - (sumX*sumX)/n);
    model[1] = -1;
    model[2] = (sumY - (model[0]*sumX))/n;
}
