#include "ransac_2Dline.h"
#include "svd.h"

int ransac_2Dline(double **data, int n, int maxT, double threshold,
					double *bestModel, int *bestInliers, int verbose) {

	if(verbose)
		printf("Start RANSAC, n=%d, maxT=%d, t=%.2lf\n\n", n, maxT, threshold);

	*bestInliers = 0;

	int T = maxT;
	int	inliers = 0;
	int Tcount = 0;
	int ndata = n;
	int nr = 2;
	int i;

	double **randSet = malloc(nr * sizeof(double *));
	if(randSet == NULL) { perror("out of memory\n"); exit(0); }
	
	for(i = 0; i < nr; i++)
	{
		randSet[i] = malloc(2 * sizeof(double));
		if(randSet[i] == NULL) { perror("out of memory\n"); exit(0); }
	}
	
	double **conSet = malloc(n * sizeof(double *));;
	if(conSet == NULL) { perror("out of memory\n"); exit(0); }
	
	for(i = 0; i < n; i++)
	{
		conSet[i] = malloc(2 * sizeof(double));
		if(conSet[i] == NULL) { perror("out of memory\n"); exit(0); }
	}
	
	double fracInliers = 0;
	double randModel[3];
	double point[2];
	double pNoOutliers = 0;
	double p = 0.99;

	srand(time(NULL)); // set rand seed

	while(T > Tcount)
	{
		if(verbose)
			printf("\n#%d ITERATION >>>>>>>>>>>>>>>>>>>>>>>\n", Tcount);

		// Select 2 points at random to form a trial model
		if(randomSelect(randSet, nr, data, &ndata) == -1)
			break;

		if(verbose)
			printf(" selected points: (%.3f, %.3f) and (%.3f, %.3f)\n", 
					randSet[0][0], randSet[0][1], randSet[1][0], 
					randSet[1][1]);

		// Fit model to the random selection of data points
		estimateModel_line(randModel, randSet, nr);

		if(verbose)
			printf(" rand model: %.3f*x + %.3f*y + %.3f = 0\n", randModel[0], 
			randModel[1], randModel[2]);

		inliers = 0;
		fracInliers = 0;

		// Evaluate distances between points and model.
		// Given a threshold, create a consensus set with the points
		// that are inliers.
		for(i = 0; i < n; i++)
		{
			point[0] = data[i][0];
			point[1] = data[i][1];

			if(fitModel_line(point, randModel, threshold))
			{
				conSet[inliers][0] = point[0];
				conSet[inliers][1] = point[1];
				inliers++;
			}
		}

		if(verbose)
			printf(" inliers = %d\n", inliers);

		if(inliers > *bestInliers)	// Largest set of inliers.
		{
			if(verbose)
				printf(" >> IT'S THE BEST MODEL !!! <<\n");
		
			// Record data for this model
			estimateModel_line(bestModel, conSet, inliers);
			*bestInliers = inliers;
			
			if(verbose)
			printf(" reestimated model: %.3f*x + %.3f*y + %.3f = 0\n",
					bestModel[0], bestModel[1], bestModel[2]);

			// Reestimate T, the number of trials to ensure we pick,
			// with probability p, a data set free of outliers.
			fracInliers = (double)inliers/n;
			pNoOutliers = 1 - pow(fracInliers, 2);
			T = log(1-p)/log(pNoOutliers);
		}

		Tcount++;
		if(Tcount > maxT)
			break;
	}

	if(bestInliers==0)
	{
		printf("\n### ERROR: ransac was unable to find a useful solution.\n");
		return(-1);
	}

	return(0);
}

int randomSelect(double **sel, int nsel, double **data, int *ndata) {

    int r = 0;
    int k = *ndata;

	if(nsel > *ndata)
	{
		printf("randomSelect: unable to select %d points from dataset[%d]\n", 
				nsel, *ndata);
		return -1;
	}

	for(int i = 0; i < nsel; i++, k--)
    {
        r = rand()%(k);

        sel[i][0] = data[r][0];
        sel[i][1] = data[r][1];

        data[r][0] = data[k-1][0];
        data[r][1] = data[k-1][1];
        
        data[k-1][0] = sel[i][0];
        data[k-1][1] = sel[i][1];
    }

	*ndata = k;
	//printf("ndata = %d\n", *ndata);

	return 0;
}

int fitModel_line(double *point, double *l, double threshold) {
    // Estimate distance between point and model
    // d = abs(a*x + b*y + c)/sqrt(a^2 + b^2)

    double d=0;

    d = fabs(l[0]*point[0] + l[1]*point[1] + l[2])/sqrt(pow(l[0], 2) + pow(l[1], 2));

    if(d<=threshold)
        return 1;
    else
        return 0;
}

void estimateModel_line(double *l, double **P, int n) {
   	int i;
   
	if(n<2) {
		perror("Need at least two points\n");
		return;
	}
	
	double p[2] = {0, 0};
	double **Q;
	double **V;
	double *W;
	Q = malloc(n * sizeof(double *));
	if(Q == NULL) { perror("out of memory\n"); exit(0); }
	
	V = malloc(n * sizeof(double *));
	if(V == NULL) { perror("out of memory\n"); exit(0); }
	
	W = malloc(2 * sizeof(double));
	if(W == NULL) { perror("out of memory\n"); exit(0); }
	
	for(i = 0; i < n; i++)
	{
		Q[i] = malloc(2 * sizeof(double));
		if(Q[i] == NULL) { perror("out of memory\n"); exit(0); }
		
		V[i] = malloc(2 * sizeof(double));
		if(V[i] == NULL) { perror("out of memory\n"); exit(0); }
	}
	
	// one = ones(m, 1);
	// centroid of all the points
	// p = (P' * one) / m;
	for(i = 0; i < n; i++)
	{
		p[0] += P[i][0];
		p[1] += P[i][1];
	}
	p[0] = p[0]/n;
	p[1] = p[1]/n;
	
	// matrix of centered coordinates
	// Q = P - one * p';
	for(i = 0; i < n; i++)
	{
		Q[i][0] = P[i][0] - p[0];
		Q[i][1] = P[i][1] - p[1];
	}
	
	// [U Sigma V] = svd(Q);
	svdcmp(Q, n, 2, W, V);
	
	// printf("U=\n"); printMatrix(Q,n,2);
	// printf("W=\n"); printVector(W,2);
	// printf("V=\n"); printMatrix(V,2,2);
	
	// the line normal is the second column of V
	// n = V(:, 2);
	// assemble the three line coefficients into a column vector
	// l = [n ; p' * n];
	l[0] = V[0][0];
	l[1] = V[0][1];
	l[2] = -(p[0]*l[0] + p[1]*l[1]);
	
	// the smallest singular value of Q
	// measures the residual fitting error
	// residue = Sigma(2, 2);
	//residue = W[1];
	
	for(i = 0; i < n; i++)
	{
		free(Q[i]);
		free(V[i]);
	}
	free(Q);
	free(V);
	free(W);
}
