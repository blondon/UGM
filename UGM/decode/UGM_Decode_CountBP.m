function [nodeLabels] = UGM_Decode_CountBP(nodePot,edgePot,edgeStruct,nodeCount,edgeCount)

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
	convTol = 1e-10;
end

% Momentum, for damping
if isfield(edgeStruct,'momentum')
	momentum = edgeStruct.momentum;
else
	momentum = .9;
end

if edgeStruct.useMex
    nodeBel = UGM_Decode_CountBPC(...
		nodePot,edgePot,nodeCount,edgeCount,...
		edgeStruct.edgeEnds,edgeStruct.nStates,edgeStruct.V,edgeStruct.E,...
		int32(edgeStruct.maxIter),momentum,convTol);
else
	nNodes = edgeStruct.nNodes;
	nEdges = edgeStruct.nEdges;
	edgeEnds = edgeStruct.edgeEnds;
	nStates = edgeStruct.nStates;
	maxStates = max(nStates);

	% Compute messages (with max-product)
	[msg_i,~] = UGM_CountBP(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,1);

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
		nodeBel(n,1:nStates(n)) = prod_of_msgs' ./ sum(prod_of_msgs);
	end
end

%% Decoding
[~,nodeLabels] = max(nodeBel,[],2);

