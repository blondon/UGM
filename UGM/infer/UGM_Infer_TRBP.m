function [nodeBel,edgeBel,logZ,H] = UGM_Infer_TRBP(nodePot,edgePot,edgeStruct,edgeDist)

if nargin >= 4
	mu = edgeDist;
elseif isfield(edgeStruct,'edgeDist')
	mu = edgeStruct.edgeDist;
else
	mu = 1;
end

% Find tree weights
if isscalar(mu)
	% Weights not provided, construct them using method specified in mu.
	mu = UGM_makeEdgeDistribution(edgeStruct,mu);
end

% Run TRBP
if edgeStruct.useMex
    [nodeBel,edgeBel,logZ,H] = UGM_Infer_TRBPC(nodePot,edgePot,edgeStruct.edgeEnds,edgeStruct.nStates,edgeStruct.V,edgeStruct.E,int32(edgeStruct.maxIter),mu);
else
	%% Non-mex version
	
	[nNodes,maxState] = size(nodePot);
	nEdges = size(edgePot,3);
	edgeEnds = edgeStruct.edgeEnds;
	nStates = double(edgeStruct.nStates);

	maximize = 0;
	new_msg = UGM_TRBP(nodePot,edgePot,edgeStruct,maximize,mu);

	% Compute nodeBel
	nodeBel = zeros(nNodes,maxState);
	for n = 1:nNodes
		% nodeBel = nodePot * (product of exponentiated messages)
		prod_of_msgs = nodePot(n,1:nStates(n))';
		for e = UGM_getEdges(n,edgeStruct)
			if n == edgeEnds(e,2)
				prod_of_msgs = prod_of_msgs .* new_msg(1:nStates(n),e).^mu(e);
			else
				prod_of_msgs = prod_of_msgs .* new_msg(1:nStates(n),e+nEdges).^mu(e);
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
	% Clamp to [0,1] (just in case)
	nodeBel(nodeBel<0) = 0; nodeBel(nodeBel>1) = 1;

	% Compute edgeBel
	if nargout > 1
		edgeBel = zeros(maxState,maxState,nEdges);
		for e = 1:nEdges
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);

			% temp1 = nodePot by all messages to n1 except from n2
			temp1 = nodePot(n1,1:nStates(n1))';
			for e2 = UGM_getEdges(n1,edgeStruct)
				if n1 == edgeEnds(e2,2)
					incoming = new_msg(1:nStates(n1),e2);
				else
					incoming = new_msg(1:nStates(n1),e2+nEdges);
				end
				if e ~= e2
					temp1 = temp1 .* incoming.^mu(e2);
				else
					temp1 = temp1 ./ incoming.^(1-mu(e2));
				end
			end

			% temp2 = nodePot by all messages to n2 except from n1
			temp2 = nodePot(n2,1:nStates(n2))';
			for e2 = UGM_getEdges(n2,edgeStruct)
				if n2 == edgeEnds(e2,2)
					incoming = new_msg(1:nStates(n2),e2);
				else
					incoming = new_msg(1:nStates(n2),e2+nEdges);
				end
				if e ~= e2
					temp2 = temp2 .* incoming.^mu(e2);
				else
					temp2 = temp2 ./ incoming.^(1-mu(e2));
				end
			end
			
			% Unnormalized edge beliefs
			eb = repmat(temp1,[1 nStates(n2)]) .* repmat(temp2',[nStates(n1) 1]) ...
				.* edgePot(1:nStates(n1),1:nStates(n2),e).^(1/mu(e));
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
	% Clamp to [0,1] (just in case)
	edgeBel(edgeBel<0) = 0; edgeBel(edgeBel>1) = 1;

	% Compute Free Energy
	if nargout > 2

		Energy1 = 0; Energy2 = 0;
		Entropy1 = 0; Entropy2 = 0;
		nodeBel = nodeBel+eps;
		edgeBel = edgeBel+eps;
		for n = 1:nNodes
			edges = UGM_getEdges(n,edgeStruct);
			nb = nodeBel(n,1:nStates(n));
			np = nodePot(n,1:nStates(n));

			% Node Entropy (note: different weighting than in Bethe)
			nodeEntropy = nb .* log(nb);
			nodeEntropy(~isfinite(nodeEntropy)) = 0;
			Entropy1 = Entropy1 - (1-sum(mu(edges)))*sum(nodeEntropy);

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

			% Pairwise Entropy (note: different weighting than in Bethe)
			edgeEntropy = eb .* log(eb);
			edgeEntropy(~isfinite(edgeEntropy)) = 0;
			Entropy2 = Entropy2 - mu(e)*sum(edgeEntropy(:));

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