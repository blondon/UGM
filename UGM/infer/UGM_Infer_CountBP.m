function [nodeBel,edgeBel,logZ,H] = UGM_Infer_CountBP(nodePot,edgePot,edgeStruct,nodeCount,edgeCount)

if nargin < 5
	if isfield(edgeStruct,'nodeCount')
		nodeCount = edgeStruct.nodeCount;
	else
		error('UGM_Infer_CountBP: requires edgeStruct.nodeCount or as 4th argument.')
	end
	if isfield(edgeStruct,'edgeCount')
		edgeCount = edgeStruct.edgeCount;
	else
		error('UGM_Infer_CountBP: requires edgeStruct.edgeCount or as 5th argument.')
	end
end

% Convergence tolerance
if isfield(edgeStruct,'convTol')
	convTol = edgeStruct.convTol;
else
	convTol = 1e-4;
end

% Momentum, for damping
if isfield(edgeStruct,'momentum')
	momentum = edgeStruct.momentum;
else
	momentum = 1;
end

if edgeStruct.useMex
    [nodeBel,edgeBel,logZ,H] = UGM_Infer_CountBPC(...
		nodePot,edgePot,nodeCount,edgeCount,...
		edgeStruct.edgeEnds,edgeStruct.nStates,edgeStruct.V,edgeStruct.E,...
		int32(edgeStruct.maxIter),momentum,convTol);
else
	%% Non-mex version

	nNodes = edgeStruct.nNodes;
	nEdges = edgeStruct.nEdges;
	edgeEnds = edgeStruct.edgeEnds;
	nStates = edgeStruct.nStates;
	maxStates = max(nStates);

	% Compute messages
	[msg_i,msg_o] = UGM_CountBP(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,0);

	% Compute nodeBel
	nodeBel = zeros(nNodes,maxStates);
	for n = 1:nNodes
		edges = UGM_getEdges(n,edgeStruct);
		prod_of_msgs = nodePot(n,1:nStates(n))';
		for e = edges
			if n == edgeEnds(e,1)
				prod_of_msgs = prod_of_msgs .* msg_i(1:nStates(n),e);
			else
				prod_of_msgs = prod_of_msgs .* msg_i(1:nStates(n),e+nEdges);
			end
		end
		
		% Safe normalize
		Z = sum(prod_of_msgs);
		if Z == 0
			nodeBel(n,1:nStates(n)) = 1 / nStates(n);
		else
			nodeBel(n,1:nStates(n)) = prod_of_msgs' ./ Z;
		end
	end
% 	if any(~isfinite(nodeBel(:)))
% 		1
% 	end
	% Clamp to [0,1] (just in case)
	nodeBel(nodeBel<0) = 0; nodeBel(nodeBel>1) = 1;

	if nargout > 1
		% Compute edge beliefs
		edgeBel = zeros(maxStates,maxStates,nEdges);
		for e = 1:nEdges
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);
			eb = edgePot(1:nStates(n1),1:nStates(n2),e).^(1/edgeCount(e)) .* ...
				 (msg_o(1:nStates(n1),e) * msg_o(1:nStates(n2),e+nEdges)');
			eb(~isfinite(eb)) = 0;
			
			% Safe normalize
			Z = sum(eb(:));
			if Z == 0
				edgeBel(1:nStates(n1),1:nStates(n2),e) = 1 / (nStates(n1)*nStates(n2));
			else
				edgeBel(1:nStates(n1),1:nStates(n2),e) = eb ./ Z;
			end
		end
	end
% 	if any(~isfinite(edgeBel(:)))
% 		1
% 	end
	% Clamp to [0,1] (just in case)
	edgeBel(edgeBel<0) = 0; edgeBel(edgeBel>1) = 1;

	if nargout > 2
		% Compute Bethe free energy
		Energy1 = 0; Energy2 = 0;
		Entropy1 = 0; Entropy2 = 0;
		nodeBel = nodeBel+eps;
		edgeBel = edgeBel+eps;
		for n = 1:nNodes
			nb = nodeBel(n,1:nStates(n));
			np = nodePot(n,1:nStates(n));
			
			% Node Entropy (can get divide by zero if beliefs at 0)
			nodeEntropy = nb .* log(nb);
			nodeEntropy(~isfinite(nodeEntropy)) = 0;
			Entropy1 = Entropy1 - nodeCount(n)*sum(nodeEntropy);

			% Node Energy
			% Note: can be infinite if (nb > 1e-10) and (np < 1e-10)
			nodeEnergy = nb .* log(np);
			nodeEnergy((nb<1e-10)&(np<1e-10)) = 0;
			Energy1 = Energy1 + sum(nodeEnergy);
		end
		for e = 1:nEdges
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);
			eb = edgeBel(1:nStates(n1),1:nStates(n2),e);
			ep = edgePot(1:nStates(n1),1:nStates(n2),e);

			% Pairwise Entropy (can get divide by zero if beliefs at 0)
			edgeEntropy = eb .* log(eb);
			edgeEntropy(~isfinite(edgeEntropy)) = 0;
			Entropy2 = Entropy2 - edgeCount(e)*sum(edgeEntropy(:));

			% Pairwise Energy
			% Note: can be infinite if (eb > 1e-10) and (ep < 1e-10)
			edgeEnergy = eb .* log(ep);
			edgeEnergy((eb<1e-10)&(ep<1e-10)) = 0;
			Energy2 = Energy2 + sum(edgeEnergy(:));
		end
		H = Entropy1 + Entropy2;
		logZ = Energy1 + Energy2 + H;
		
	end

end
