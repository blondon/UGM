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
	
	double *nodePot,*edgePot,*nodeCount,*edgeCount,
			*nodeBel,
			q1,q2,d1,d2,
			*powEdgePot,
			*msg_i,*msg_o,*tmp_i,*tmp_o,*old_msg_i,*old_msg_o,
			z,momentum,sumAbsDiff,convTol;
	
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
	momentum = ((double*)mxGetPr(prhs[9]))[0];
	convTol = ((double*)mxGetPr(prhs[10]))[0];
	
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
	
	/* Local copy of edgePot^(1/edgeCount) */
	powEdgePot = mxCalloc(maxState*maxState*nEdges,sizeof(double));
	for (e = 0; e < nEdges; e++)
	{
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		for (s1 = 0; s1 < nStates[n1]; s1++)
			for (s2 = 0; s2 < nStates[n1]; s2++)
				powEdgePot[s1+maxState*(s2+maxState*e)] = pow(edgePot[s1+maxState*(s2+maxState*e)], 1.0/edgeCount[e]);
	}
	
	/* Loop variables */
	nMessages = maxState * nEdges * 2;
	msg_i = mxCalloc(nMessages,sizeof(double));
	msg_o = mxCalloc(nMessages,sizeof(double));
	tmp_i = mxCalloc(nMessages,sizeof(double));
	tmp_o = mxCalloc(nMessages,sizeof(double));
	old_msg_i = mxCalloc(nMessages,sizeof(double));
	old_msg_o = mxCalloc(nMessages,sizeof(double));

	/* Initialize */
	for (e = 0; e < nEdges; e++)
	{
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		/* Init messages */
		for (s = 0; s < nStates[n1]; s++) {
			msg_i[s+maxState*e] = 1.0/nStates[n1];
			msg_o[s+maxState*e] = 1.0/nStates[n1];
			old_msg_i[s+maxState*e] = 1.0/nStates[n1];
			old_msg_o[s+maxState*e] = 1.0/nStates[n1];
		}
		for (s = 0; s < nStates[n2]; s++) {
			msg_i[s+maxState*(e+nEdges)] = 1.0/nStates[n2];
			msg_o[s+maxState*(e+nEdges)] = 1.0/nStates[n2];
			old_msg_i[s+maxState*(e+nEdges)] = 1.0/nStates[n2];
			old_msg_o[s+maxState*(e+nEdges)] = 1.0/nStates[n2];
		}
	}
	
	/* Main loop */
	for (iter = 0; iter < maxIter; iter++)
	{
		/* Temp messages */
		for (e = 0; e < nEdges; e++)
		{
			n1 = edgeEnds[e]-1;
			n2 = edgeEnds[e+nEdges]-1;
			
			/* Incoming */
			for (s1 = 0; s1 < nStates[n1]; s1++) {
				tmp_i[s1+maxState*e] = 0.0;
				for (s2 = 0; s2 < nStates[n2]; s2++) {
					z = powEdgePot[s1+maxState*(s2+maxState*e)] * msg_o[s2+maxState*(e+nEdges)];
					if (tmp_i[s1+maxState*e] < z)
						tmp_i[s1+maxState*e] = z;
				}
			}
			for (s2 = 0; s2 < nStates[n2]; s2++) {
				tmp_i[s2+maxState*(e+nEdges)] = 0.0;
				for (s1 = 0; s1 < nStates[n1]; s1++) {
					z = powEdgePot[s1+maxState*(s2+maxState*e)] * msg_o[s1+maxState*e];
					if (tmp_i[s2+maxState*(e+nEdges)] < z)
						tmp_i[s2+maxState*(e+nEdges)] = z;
				}
			}

			/* Outgoing */
			for (s = 0; s < nStates[n1]; s++)
				tmp_o[s+maxState*e] = nodePot[n1+nNodes*s];
			for(Vind = V[n1]-1; Vind < V[n1+1]-1; Vind++)
			{
				e2 = E[Vind] - 1;
				if (e != e2)
				{
					if (n1 == edgeEnds[e2]-1) {
						for (s = 0; s < nStates[n1]; s++)
							tmp_o[s+maxState*e] *= msg_i[s+maxState*e2];
					}
					else {
						for (s = 0; s < nStates[n1]; s++)
							tmp_o[s+maxState*e] *= msg_i[s+maxState*(e2+nEdges)];
					}
				}
			}
			for (s = 0; s < nStates[n2]; s++)
				tmp_o[s+maxState*(e+nEdges)] = nodePot[n2+nNodes*s];
			for(Vind = V[n2]-1; Vind < V[n2+1]-1; Vind++)
			{
				e2 = E[Vind] - 1;
				if (e != e2)
				{
					if (n2 == edgeEnds[e2]-1) {
						for (s = 0; s < nStates[n2]; s++)
							tmp_o[s+maxState*(e+nEdges)] *= msg_i[s+maxState*e2];
					}
					else {
						for (s = 0; s < nStates[n2]; s++)
							tmp_o[s+maxState*(e+nEdges)] *= msg_i[s+maxState*(e2+nEdges)];
					}
				}
			}
		}
		
		/* New messages */
		for (e = 0; e < nEdges; e++)
		{
			n1 = edgeEnds[e]-1;
			n2 = edgeEnds[e+nEdges]-1;
			
			/* Exponent variables */
			q1 = (1.0-nodeCount[n1]) / (V[n1+1]-V[n1]);
			q2 = (1.0-nodeCount[n2]) / (V[n2+1]-V[n2]);
			d1 = edgeCount[e] - q1 + 1.0;
			d2 = edgeCount[e] - q2 + 1.0;

			/* Incoming */
			z = 0.0;
			for (s = 0; s < nStates[n1]; s++) {
				idx = s+maxState*e;
				msg_i[idx] = zeroIfNaN(
								pow(tmp_i[idx], edgeCount[e]/d1)
							  * pow(tmp_o[idx], (q1-edgeCount[e])/d1)
							 );
				z += msg_i[idx];
			}
			if (z != 0.0) /* Normalize */
			{
				for (s = 0; s < nStates[n1]; s++)
					msg_i[s+maxState*e] = msg_i[s+maxState*e] / z;
			}
			/*
			else
			{
				for (s = 0; s < nStates[n1]; s++)
					msg_i[s+maxState*e] = 1.0 / nStates[n1];
			}
			*/
			
			z = 0.0;
			for (s = 0; s < nStates[n2]; s++) {
				idx = s+maxState*(e+nEdges);
				msg_i[idx] = zeroIfNaN(
								pow(tmp_i[idx], edgeCount[e]/d2)
							  * pow(tmp_o[idx], (q2-edgeCount[e])/d2)
							 );
				z += msg_i[idx];
			}
			if (z != 0.0) /* Normalize */
			{
				for (s = 0; s < nStates[n2]; s++)
					msg_i[s+maxState*(e+nEdges)] = msg_i[s+maxState*(e+nEdges)] / z;
			}
			/*
			else
			{
				for (s = 0; s < nStates[n2]; s++)
					msg_i[s+maxState*(e+nEdges)] = 1.0 / nStates[n2];
			}
			*/
			
			/* Outgoing */
			z = 0.0;
			for (s = 0; s < nStates[n1]; s++) {
				idx = s+maxState*e;
				msg_o[idx] = zeroIfNaN(
								pow(tmp_i[idx], (q1-1)/d1)
							  * pow(tmp_o[idx], 1/d1)
							 );
				z += msg_o[idx];
			}
			if (z != 0.0) /* Normalize */
			{
				for (s = 0; s < nStates[n1]; s++)
					msg_o[s+maxState*e] = msg_o[s+maxState*e] / z;
			}
			/*
			else
			{
				for (s = 0; s < nStates[n1]; s++)
					msg_o[s+maxState*e] = 1.0 / nStates[n1];
			}
			*/
			
			z = 0.0;
			for (s = 0; s < nStates[n2]; s++) {
				idx = s+maxState*(e+nEdges);
				msg_o[idx] = zeroIfNaN(
								pow(tmp_i[idx], (q2-1)/d2)
							  * pow(tmp_o[idx], 1/d2)
							 );
				z += msg_o[idx];
			}
			if (z != 0.0) /* Normalize */
			{
				for (s = 0; s < nStates[n2]; s++)
					msg_o[s+maxState*(e+nEdges)] = msg_o[s+maxState*(e+nEdges)] / z;
			}
			/*
			else
			{
				for (s = 0; s < nStates[n2]; s++)
					msg_o[s+maxState*(e+nEdges)] = 1.0 / nStates[n2];
			}
			*/
			
			/* Damping */
			if (momentum < 1.0)
			{
				for (s = 0; s < nStates[n1]; s++) {
					idx = s+maxState*e;
					msg_i[idx] = pow(old_msg_i[idx],1-momentum) * pow(msg_i[idx],momentum);
					msg_o[idx] = pow(old_msg_o[idx],1-momentum) * pow(msg_o[idx],momentum);
				}
				for (s = 0; s < nStates[n2]; s++) {
					idx = s+maxState*(e+nEdges);
					msg_i[idx] = pow(old_msg_i[idx],1-momentum) * pow(msg_i[idx],momentum);
					msg_o[idx] = pow(old_msg_o[idx],1-momentum) * pow(msg_o[idx],momentum);
				}
			}
		}
	
		/* Check convergence */
		notConverged = 0;
		for (s = 0; s < maxState; s++) {
			for (e = 0; e < nEdges*2; e++) {
				idx = s+maxState*e;
				if (absDif(msg_i[idx],old_msg_i[idx]) >= convTol)
					notConverged++;
				if (absDif(msg_o[idx],old_msg_o[idx]) >= convTol)
					notConverged++;
				old_msg_i[idx] = msg_i[idx];
				old_msg_o[idx] = msg_o[idx];
			}
		}
		if (notConverged == 0) {
			iter++;
			break;
		}
		/* Sum-abs-diff convergence criteria is not consistent with non-mex version.
		sumAbsDiff = 0.0;
		for (idx = 0; idx < nMessages; idx++) {
			sumAbsDiff += absDif(msg_i[idx],old_msg_i[idx]);
			sumAbsDiff += absDif(msg_o[idx],old_msg_o[idx]);
			old_msg_i[idx] = msg_i[idx];
			old_msg_o[idx] = msg_o[idx];
		}
		if (sumAbsDiff < convTol) {
			iter++;
			break;
		}
		*/
	}
	
	/*
	if(iter == maxIter)
		printf("CountBP did not converge after %d iterations\n",maxIter);
	printf("Stopped after %d iterations\n",iter);
	*/
	
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
					nodeBel[n+nNodes*s] *= msg_i[s+maxState*e];
			}
			else {
				for(s = 0; s < nStates[n]; s++)
					nodeBel[n+nNodes*s] *= msg_i[s+maxState*(e+nEdges)];
			}
		}
		/* Normalize */
		z = 0.0;
		for(s = 0; s < nStates[n]; s++)
			z += nodeBel[n+nNodes*s];
		for(s = 0; s < nStates[n]; s++) {
			/* Note: unsafe normalize! */
			nodeBel[n+nNodes*s] /= z;
			/* Clamp to [0,1] (just in case) */
			if (nodeBel[n+nNodes*s] < 0.0)
				nodeBel[n+nNodes*s] = 0.0;
			else if (nodeBel[n+nNodes*s] > 1.0)
				nodeBel[n+nNodes*s] = 1.0;
		}
	}
	
	/* Free memory */
	mxFree(powEdgePot);
	mxFree(msg_i);
	mxFree(msg_o);
	mxFree(tmp_i);
	mxFree(tmp_o);
	mxFree(old_msg_i);
	mxFree(old_msg_o);
}
