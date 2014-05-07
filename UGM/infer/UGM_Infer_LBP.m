function [nodeBel,edgeBel,logZ,H] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct)

if edgeStruct.useMex
    [nodeBel,edgeBel,logZ,H] = UGM_Infer_LBPC(nodePot,edgePot,edgeStruct.edgeEnds,edgeStruct.nStates,edgeStruct.V,edgeStruct.E,int32(edgeStruct.maxIter));
else
	%% Non-mex version
	
	[nNodes,maxState] = size(nodePot);
	nEdges = size(edgePot,3);
	edgeEnds = edgeStruct.edgeEnds;
	V = edgeStruct.V;
	E = edgeStruct.E;
	nStates = edgeStruct.nStates;

	maximize = 0;
	new_msg = UGM_LoopyBP(nodePot,edgePot,edgeStruct,maximize);


	% Compute nodeBel
	nodeBel = zeros(nNodes,maxState);
	for n = 1:nNodes
		edges = UGM_getEdges(n,edgeStruct);
		prod_of_msgs = nodePot(n,1:nStates(n))';
		for e = edges
			if n == edgeEnds(e,2)
				prod_of_msgs = prod_of_msgs .* new_msg(1:nStates(n),e);
			else
				prod_of_msgs = prod_of_msgs .* new_msg(1:nStates(n),e+nEdges);
			end
		end
		nodeBel(n,1:nStates(n)) = prod_of_msgs'./sum(prod_of_msgs);
	end

	if nargout > 1
		% Compute edge beliefs
		edgeBel = zeros(maxState,maxState,nEdges);
		for e = 1:nEdges
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);
			belN1 = nodeBel(n1,1:nStates(n1))'./new_msg(1:nStates(n1),e+nEdges);
			belN2 = nodeBel(n2,1:nStates(n2))'./new_msg(1:nStates(n2),e);
			b1=repmat(belN1,1,nStates(n2));
			b2=repmat(belN2',nStates(n1),1);
			eb = b1.*b2.*edgePot(1:nStates(n1),1:nStates(n2),e);
			edgeBel(1:nStates(n1),1:nStates(n2),e) = eb./sum(eb(:));
		end
	end

	if nargout > 2
		% Compute Bethe free energy
		Energy1 = 0; Energy2 = 0;
		Entropy1 = 0; Entropy2 = 0;
		nodeBel = nodeBel+eps;
		edgeBel = edgeBel+eps;
		for n = 1:nNodes
			edges = E(V(n):V(n+1)-1);
			nNbrs = length(edges);

			% Node Entropy (can get divide by zero if beliefs at 0)
			nb = nodeBel(n,1:nStates(n));
			Entropy1 = Entropy1 - (1-nNbrs)*sum(nb.*log(nb));

			% Node Energy
			Energy1 = Energy1 + sum(nb.*log(nodePot(n,1:nStates(n))));
		end
		for e = 1:nEdges
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);

			% Pairwise Entropy (can get divide by zero if beliefs at 0)
			eb = edgeBel(1:nStates(n1),1:nStates(n2),e);
			Entropy2 = Entropy2 - sum(eb(:).*log(eb(:)));

			% Pairwise Energy
			ep = edgePot(1:nStates(n1),1:nStates(n2),e);
			Energy2 = Energy2 + sum(eb(:).*log(ep(:)));
		end
		H = Entropy1 + Entropy2;
		logZ = Energy1 + Energy2 + H;
	end
	
end