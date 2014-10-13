#include <math.h>
#include "mex.h"
#include "UGM_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Variables */
	int n,s,e,e2,n1,n2,neigh,Vind,Vind2,s1,s2,idx,
		nNodes,nEdges,maxState,nMessages,
		iter,maxIter,notConverged,nNbrs,
		*edgeEnds,*nStates,*V,*E,*y;
	
	double *nodePot,*edgePot,*logNodePot,*logEdgePot,
			*nodeCount,*edgeCount,*auxNodeCount,
			*nodeBel,
			*imsg,*omsg,*imsg_old,*omsg_old,*tmp1,*tmp2,
			z,sumAbsDiff,convTol;
		
	/* Input */
	nodePot = mxGetPr(prhs[0]);
	edgePot = mxGetPr(prhs[1]);
	nodeCount = mxGetPr(prhs[2]);
	edgeCount = mxGetPr(prhs[3]);
	edgeEnds = (int*)mxGetPr(prhs[4]);
	nStates = (int*)mxGetPr(prhs[5]);
	V = (int*)mxGetPr(prhs[6]);
	E = (int*)mxGetPr(prhs[7]);
	maxIter = ((int*)mxGetPr(prhs[8]))[0];
	convTol = ((double*)mxGetPr(prhs[9]))[0];
	
	if (!mxIsClass(prhs[4],"int32")
	||!mxIsClass(prhs[5],"int32")
	||!mxIsClass(prhs[6],"int32")
	||!mxIsClass(prhs[7],"int32")
	||!mxIsClass(prhs[8],"int32"))
		mexErrMsgTxt("edgeEnds, nStates, V, E, maxIter must be int32");
	
	/* Compute sizes */
	nNodes = mxGetDimensions(prhs[0])[0];
	maxState = mxGetDimensions(prhs[0])[1];
	nEdges = mxGetDimensions(prhs[4])[0];
	
	/* Output */
	plhs[0] = mxCreateDoubleMatrix(nNodes,maxState,mxREAL);
	nodeBel = mxGetPr(plhs[0]);

	/* Precompute log(nodePot), log(edgePot) */
	logNodePot = mxCalloc(nNodes*maxState,sizeof(double));
	for (n = 0; n < nNodes; n++) {
		for (s = 0; s < nStates[n]; s++)
			logNodePot[n+nNodes*s] = log(nodePot[n+nNodes*s]);
	}
	logEdgePot = mxCalloc(maxState*maxState*nEdges,sizeof(double));
	for (e = 0; e < nEdges; e++) {
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		for (s1 = 0; s1 < nStates[n1]; s1++) {
			for (s2 = 0; s2 < nStates[n2]; s2++)
				logEdgePot[s1+maxState*(s2+maxState*e)] = log(edgePot[s1+maxState*(s2+maxState*e)]);
		}
	}
	
	/* Precompute aux node counts */
	auxNodeCount = mxCalloc(nNodes,sizeof(double));
	for (n = 0; n < nNodes; n++) {
		auxNodeCount[n] = nodeCount[n];
		for (Vind = V[n]-1; Vind < V[n+1]-1; Vind++) {
			e = E[Vind]-1;
			auxNodeCount[n] += edgeCount[e];
		}
	}
	
	/* Loop variables */
	nMessages = maxState * nEdges * 2;
	imsg = mxCalloc(nMessages,sizeof(double));
	omsg = mxCalloc(nMessages,sizeof(double));
	imsg_old = mxCalloc(nMessages,sizeof(double));
	omsg_old = mxCalloc(nMessages,sizeof(double));
	tmp1 = mxCalloc(maxState,sizeof(double));
	tmp2 = mxCalloc(maxState,sizeof(double));

	/* Initialize */
	for (e = 0; e < nEdges; e++)
	{
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		/* Init messages */
		for (s = 0; s < nStates[n1]; s++) {
			imsg[s+maxState*e] = 0.0;
			omsg[s+maxState*e] = 0.0;
			imsg_old[s+maxState*e] = 0.0;
			omsg_old[s+maxState*e] = 0.0;
		}
		for (s = 0; s < nStates[n2]; s++) {
			imsg[s+maxState*(e+nEdges)] = 0.0;
			omsg[s+maxState*(e+nEdges)] = 0.0;
			imsg_old[s+maxState*(e+nEdges)] = 0.0;
			omsg_old[s+maxState*(e+nEdges)] = 0.0;
		}
	}
	
	
	/* Main loop */
	for (iter = 0; iter < maxIter; iter++)
	{
		/* Iterate over nodes */
		for (n = 0; n < nNodes; n++)
		{
			/* Incoming messages */
			for (Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
			{
				e = E[Vind] - 1;
				n1 = edgeEnds[e]-1;
				n2 = edgeEnds[e+nEdges]-1;
				if (n == n1)
				{
					for (s1 = 0; s1 < nStates[n1]; s1++) {
						tmp1[s1] = -INFINITY;
						for (s2 = 0; s2 < nStates[n2]; s2++) {
							z = (
									logEdgePot[s1+maxState*(s2+maxState*e)]
								  + omsg[s2+maxState*(e+nEdges)]
								) / edgeCount[e];
							if (tmp1[s1] < z)
								tmp1[s1] = z;
						}
						tmp1[s1] *= edgeCount[e];
					}
					logNormalize(tmp1,nStates[n1]);
					for (s1 = 0; s1 < nStates[n1]; s1++)
						imsg[s1+maxState*e] = tmp1[s1];
				}
				else
				{
					for (s2 = 0; s2 < nStates[n2]; s2++) {
						tmp2[s2] = -INFINITY;
						for (s1 = 0; s1 < nStates[n1]; s1++) {
							z = (
									logEdgePot[s1+maxState*(s2+maxState*e)]
								  + omsg[s1+maxState*e]
								) / edgeCount[e];
							if (tmp2[s2] < z)
								tmp2[s2] = z;
						}
						tmp2[s2] *= edgeCount[e];
					}
					logNormalize(tmp2,nStates[n2]);
					for (s2 = 0; s2 < nStates[n2]; s2++)
						imsg[s2+maxState*(e+nEdges)] = tmp2[s2];
				}
			}
			
			/* Outgoing messages */
			for (s = 0; s < nStates[n]; s++)
			{
				tmp1[s] = logNodePot[n+nNodes*s];
				for (Vind = V[n]-1; Vind < V[n+1]-1; Vind++) {
					e = E[Vind] - 1;
					if (n == edgeEnds[e]-1)
						tmp1[s] += imsg[s+maxState*e];
					else
						tmp1[s] += imsg[s+maxState*(e+nEdges)];
				}
			}
			for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
			{
				e = E[Vind] - 1;
				if (n == edgeEnds[e]-1) {
					for (s = 0; s < nStates[n]; s++)
						tmp2[s] = tmp1[s] * (edgeCount[e]/auxNodeCount[n]) - imsg[s+maxState*e];
					logNormalize(tmp2,nStates[n]);
					for (s = 0; s < nStates[n]; s++)
						omsg[s+maxState*e] = tmp2[s];
				}
				else {
					for (s = 0; s < nStates[n]; s++)
						tmp2[s] = tmp1[s] * (edgeCount[e]/auxNodeCount[n]) - imsg[s+maxState*(e+nEdges)];
					logNormalize(tmp2,nStates[n]);
					for (s = 0; s < nStates[n]; s++)
						omsg[s+maxState*(e+nEdges)] = tmp2[s];
				}
			}
		}
	
		/* Check convergence */
		notConverged = 0;
		for (s = 0; s < maxState; s++) {
			for (e = 0; e < nEdges*2; e++) {
				idx = s+maxState*e;
				if (absDif(imsg[idx],imsg_old[idx]) >= convTol)
					notConverged = 1;
				if (absDif(omsg[idx],omsg_old[idx]) >= convTol)
					notConverged = 1;
				imsg_old[idx] = imsg[idx];
				omsg_old[idx] = omsg[idx];
			}
		}
		if (notConverged == 0) {
			iter++;
			break;
		}
	}
	
	/*
	if(iter == maxIter)
		printf("CountBP did not converge after %d iterations\n",maxIter);
	printf("Stopped after %d iterations\n",iter);
	*/
	
	/* Convert to exponential space */
	for (e = 0; e < nEdges; e++)
	{
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		for (s = 0; s < nStates[n1]; s++) {
			imsg[s+maxState*e] = exp(imsg[s+maxState*e]);
			omsg[s+maxState*e] = exp(omsg[s+maxState*e]);
		}
		for (s = 0; s < nStates[n2]; s++) {
			imsg[s+maxState*(e+nEdges)] = exp(imsg[s+maxState*(e+nEdges)]);
			omsg[s+maxState*(e+nEdges)] = exp(omsg[s+maxState*(e+nEdges)]);
		}
	}
	
	/* Compute nodeBel */
	for (n = 0; n < nNodes; n++)
	{
		/* Init to nodePot */
		for (s = 0; s < nStates[n]; s++)
			nodeBel[n+nNodes*s] = nodePot[n+nNodes*s];
		/* Multiply by product of incoming messages */
		for (Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
		{
			e = E[Vind]-1;
			n1 = edgeEnds[e]-1;
			n2 = edgeEnds[e+nEdges]-1;
			if (n == n1) {
				for(s = 0; s < nStates[n]; s++)
					nodeBel[n+nNodes*s] *= imsg[s+maxState*e];
			}
			else {
				for(s = 0; s < nStates[n]; s++)
					nodeBel[n+nNodes*s] *= imsg[s+maxState*(e+nEdges)];
			}
		}
		/* Normalize */
		z = 0.0;
		for(s = 0; s < nStates[n]; s++)
			z += nodeBel[n+nNodes*s];
		for(s = 0; s < nStates[n]; s++) {
			if (z == 0.0) {
				/* Uniform beliefs */
				nodeBel[n+nNodes*s] = 1.0 / nStates[n];
			}
			else {
				/* Safe to normalize */
				nodeBel[n+nNodes*s] /= z;
				/* Clamp to [0,1] (just in case) */
				if (nodeBel[n+nNodes*s] < 0.0)
					nodeBel[n+nNodes*s] = 0.0;
				else if (nodeBel[n+nNodes*s] > 1.0)
					nodeBel[n+nNodes*s] = 1.0;
			}
		}
	}
	
	/* Free memory */
	mxFree(logNodePot);
	mxFree(logEdgePot);
	mxFree(imsg);
	mxFree(omsg);
	mxFree(imsg_old);
	mxFree(omsg_old);
	mxFree(tmp1);
	mxFree(tmp2);
	mxFree(auxNodeCount);
}
