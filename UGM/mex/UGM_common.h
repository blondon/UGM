
int getMaxState(int *nStates,int nNodes)
{
   int n, maxState=0;
   
   for(n = 0; n < nNodes; n++)
   {
      if(nStates[n] > maxState)
         maxState = nStates[n];
   }
   return maxState;
}

double absDif(double a, double b)
 {
     if (a > b)
     {
         return a-b;
     }
     else
     {
         return b-a;
     }
 }
 
double zeroIfNaN(double x)
{
	if (mxIsFinite(x))
	{
		return x;
	}
	return 0.0;
}

void logNormalize(double *x, int length)
{
	double maxval = -INFINITY;
	for (int i = 0; i < length; i++) {
		if (maxval < x[i])
			maxval = x[i];
	}
	for (int i = 0; i < length; i++) {
		x[i] -= maxval;
	}
}

double logSumExp(double *x, int length)
{
	double maxval = -INFINITY;
	for (int i = 0; i < length; i++) {
		if (maxval < x[i])
			maxval = x[i];
	}
	double sum = 0.0;
	for (int i = 0; i < length; i++) {
		sum += exp(x[i] - maxval);
	}
	return log(sum) + maxval;
}

 
