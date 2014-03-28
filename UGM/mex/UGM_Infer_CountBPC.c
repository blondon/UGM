#include <math.h>
#include "mex.h"
#include "UGM_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Variables */
	int n,s,e,e2,n1,n2,neigh,Vind,Vind2,s1,s2,
			nNodes,nEdges,maxState,dims[3],
			iter,maxIter,nNbrs,notConverged,
			*edgeEnds,*nStates,*V,*E,*y;
	
	double *nodePot,*edgePot,*nodeCount,*edgeCount,
			*nodeBel,*edgeBel,*logZ,
			q1,q2,d1,d2,
			*msg_i,*msg_o,*tmp_i,*tmp_o,*old_msg_i,*old_msg_o,
			z,convTol,
			energy1,energy2,entropy1,entropy2;
// 			*prodMsgs,*oldMsgs,*newMsgs,*tmp;
	
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
	dims[0] = maxState;
	dims[1] = maxState;
	dims[2] = nEdges;
	plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	nodeBel = mxGetPr(plhs[0]);
	edgeBel = mxGetPr(plhs[1]);
	logZ = mxGetPr(plhs[2]);
	
	/* Loop variables */
	msg_i = mxCalloc(maxState*nEdges*2,sizeof(double));
	msg_o = mxCalloc(maxState*nEdges*2,sizeof(double));
	tmp_i = mxCalloc(maxState*nEdges*2,sizeof(double));
	tmp_o = mxCalloc(maxState*nEdges*2,sizeof(double));
	old_msg_i = mxCalloc(maxState*nEdges*2,sizeof(double));
	old_msg_o = mxCalloc(maxState*nEdges*2,sizeof(double));
	convTol = 1e-4;

	/* Initialize */
	for (e = 0; e < nEdges; e++)
	{
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		/* Init messages */
		for (s = 0; s < nStates[n1]; s++) {
			msg_i[s+maxState*e] = 1.0/nStates[n1];
			msg_o[s+maxState*e] = 1.0/nStates[n1];
		}
		for (s = 0; s < nStates[n2]; s++) {
			msg_i[s+maxState*(e+nEdges)] = 1.0/nStates[n2];
			msg_o[s+maxState*(e+nEdges)] = 1.0/nStates[n2];
		}
		/* Raise edgePot to power (1/edgeCount[e]) */
		for (s1 = 0; s1 < nStates[n1]; s1++) {
			for (s2 = 0; s2 < nStates[n1]; s2++)
				edgePot[s1+maxState*(s2+maxState*e)] = pow(edgePot[s1+maxState*(s2+maxState*e)], 1.0/edgeCount[e]);
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
				for (s2 = 0; s2 < nStates[n1]; s2++)
					tmp_i[s1+maxState*e] += edgePot[s1+maxState*(s2+maxState*e)] * msg_o[s2+maxState*(e+nEdges)];
			}
			for (s2 = 0; s2 < nStates[n1]; s2++) {
				tmp_i[s2+maxState*(e+nEdges)] = 0.0;
				for (s1 = 0; s1 < nStates[n1]; s1++)
					tmp_i[s2+maxState*(e+nEdges)] += edgePot[s1+maxState*(s2+maxState*e)] * msg_o[s1+maxState*e];
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
			q1 = (1-nodeCount[n1]) / (V[n1]-V[n1+1]);
			q2 = (1-nodeCount[n2]) / (V[n2]-V[n2+1]);
			d1 = edgeCount[e] - q1 + 1;
			d2 = edgeCount[e] - q2 + 1;

			/* Incoming */
			z = 0.0;
			for (s = 0; s < nStates[n1]; s++) {
				msg_i[s+maxState*e] = pow(tmp_i[s+maxState*e],edgeCount[e]/d1)
									* pow(tmp_o[s+maxState*e],(q1-edgeCount[e])/d1);
				z += msg_i[s+maxState*e];
			}
			for (s = 0; s < nStates[n1]; s++)
				msg_i[s+maxState*e] /= z;
			z = 0.0;
			for (s = 0; s < nStates[n2]; s++) {
				msg_i[s+maxState*(e+nEdges)] = pow(tmp_i[s+maxState*(e+nEdges)],edgeCount[e]/d2)
											* pow(tmp_o[s+maxState*(e+nEdges)],(q2-edgeCount[e])/d2);
				z += msg_i[s+maxState*(e+nEdges)];
			}
			for (s = 0; s < nStates[n2]; s++)
				msg_i[s+maxState*(e+nEdges)] /= z;
			
			/* Outgoing */
			z = 0.0;
			for (s = 0; s < nStates[n1]; s++) {
				msg_o[s+maxState*e] = pow(tmp_i[s+maxState*e],(q1-1)/d1)
									* pow(tmp_o[s+maxState*e],1/d1);
				z += msg_o[s+maxState*e];
			}
			for (s = 0; s < nStates[n1]; s++)
				msg_o[s+maxState*e] /= z;
			z = 0.0;
			for (s = 0; s < nStates[n2]; s++) {
				msg_o[s+maxState*(e+nEdges)] = pow(tmp_i[s+maxState*(e+nEdges)],(q2-1)/d2)
											* pow(tmp_o[s+maxState*(e+nEdges)],1/d2);
				z += msg_o[s+maxState*(e+nEdges)];
			}
			for (s = 0; s < nStates[n2]; s++)
				msg_o[s+maxState*(e+nEdges)] /= z;
		}
		
		/* Check convergence */
		notConverged = 0;
		for (s = 0; s < maxState; s++) {
			for (e = 0; e < nEdges*2; e++) {
				if (absDif(msg_i[s+maxState*e],old_msg_i[s+maxState*e]) >= convTol)
					notConverged++;
				if (absDif(msg_o[s+maxState*e],old_msg_o[s+maxState*e]) >= convTol)
					notConverged++;
				old_msg_i[s+maxState*e] = msg_i[s+maxState*e];
				old_msg_o[s+maxState*e] = msg_o[s+maxState*e];
			}
		}
		if (notConverged == 0) {
			iter++;
			break;
		}
	}
	
// 	if(iter == maxIter)
// 		printf("LBP reached maxIter of %d iterations\n",maxIter);
// 	printf("Stopped after %d iterations\n",iter);
	

// 	/* Compute nodeBel */
// 	for(n = 0; n < nNodes; n++)
// 	{
// 		for(s = 0; s < nStates[n]; s++)
// 			prodMsgs[s+maxState*n] = nodePot[n+nNodes*s];
// 		
// 		for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
// 		{
// 			e = E[Vind]-1;
// 			n1 = edgeEnds[e]-1;
// 			n2 = edgeEnds[e+nEdges]-1;
// 			
// 			if (n == n2)
// 			{
// 				for(s = 0; s < nStates[n]; s++)
// 				{
// 					prodMsgs[s+maxState*n] *= newMsgs[s+maxState*e];
// 				}
// 			}
// 			else
// 			{
// 				for(s = 0; s < nStates[n]; s++)
// 				{
// 					prodMsgs[s+maxState*n] *= newMsgs[s+maxState*(e+nEdges)];
// 				}
// 			}
// 		}
// 		
// 		z = 0;
// 		for(s = 0; s < nStates[n]; s++)
// 		{
// 			nodeBel[n + nNodes*s] = prodMsgs[s+maxState*n];
// 			z = z + nodeBel[n+nNodes*s];
// 		}
// 		for(s = 0; s < nStates[n]; s++)
// 			nodeBel[n + nNodes*s] /= z;
// 	}
// 	
// 	
// 	/* Compute edgeBel */
// 	for(e = 0; e < nEdges; e++)
// 	{
// 		n1 = edgeEnds[e]-1;
// 		n2 = edgeEnds[e+nEdges]-1;
// 		z = 0;
// 		for(s1 = 0; s1 < nStates[n1]; s1++)
// 		{
// 			for(s2 = 0; s2 < nStates[n2]; s2++)
// 			{
// 				edgeBel[s1+maxState*(s2+maxState*e)] = nodeBel[n1+nNodes*s1]/newMsgs[s1+maxState*(e+nEdges)];
// 				edgeBel[s1+maxState*(s2+maxState*e)] *= nodeBel[n2+nNodes*s2]/newMsgs[s2+maxState*e];
// 				edgeBel[s1+maxState*(s2+maxState*e)] *= edgePot[s1+maxState*(s2+maxState*e)];
// 				z += edgeBel[s1+maxState*(s2+maxState*e)];
// 			}
// 		}
// 		for(s1 = 0; s1 < nStates[n1]; s1++)
// 		{
// 			for(s2 = 0; s2 < nStates[n2]; s2++)
// 				edgeBel[s1+maxState*(s2+maxState*e)] /= z;
// 		}
// 	}
// 	
// 	/* Compute Bethe Free Energy */
// 	energy1 = 0;
// 	energy2 = 0;
// 	entropy1 = 0;
// 	entropy2 = 0;
// 	for(n = 0; n < nNodes; n++)
// 	{
// 		nNbrs = V[n+1]-V[n];
// 		for(s = 0; s < nStates[n]; s++)
// 		{
// 			if(nodeBel[n+nNodes*s] > 1e-10)
// 				entropy1 += (nNbrs-1)*nodeBel[n+nNodes*s]*log(nodeBel[n+nNodes*s]);
// 			
// 			energy1 -= nodeBel[n+nNodes*s]*log(nodePot[n+nNodes*s]);
// 		}
// 	}
// 	for(e = 0; e < nEdges; e++)
// 	{
// 		n1 = edgeEnds[e]-1;
// 		n2 = edgeEnds[e+nEdges]-1;
// 		
// 		for(s1 = 0; s1 < nStates[n1];s1++)
// 		{
// 			for(s2 = 0; s2 < nStates[n2]; s2++)
// 			{
// 				if(edgeBel[s1+maxState*(s2+maxState*e)] > 1e-10)
// 				{
// 					entropy2 -= edgeBel[s1+maxState*(s2+maxState*e)]*log(edgeBel[s1+maxState*(s2+maxState*e)]);
// 				}
// 				energy2 -= edgeBel[s1+maxState*(s2+maxState*e)]*log(edgePot[s1+maxState*(s2+maxState*e)]);
// 			}
// 		}
// 	}
// 	logZ[0] = -energy1-energy2+entropy1+entropy2;
	
	
	/* Free memory */
	mxFree(msg_i);
	mxFree(msg_o);
	mxFree(tmp_i);
	mxFree(tmp_o);
	mxFree(old_msg_i);
	mxFree(old_msg_o);
// 	mxFree(prodMsgs);
// 	mxFree(oldMsgs);
// 	mxFree(newMsgs);
// 	mxFree(tmp);
}
