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
	%% Non-mex version

	nNodes = double(edgeStruct.nNodes);
	nEdges = double(edgeStruct.nEdges);
	edgeEnds = double(edgeStruct.edgeEnds);
	nStates = double(edgeStruct.nStates);
	maxState = max(nStates);
    
	% Compute messages
	[imsg,~] = UGM_CountBP(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,convTol,1);

    % Compute nodeBel
	nodeBel = zeros(nNodes,maxState);
	for n = 1:nNodes
		prod_of_msgs = nodePot(n,1:nStates(n))';
		for e = UGM_getEdges(n,edgeStruct);
			if n == edgeEnds(e,1)
				prod_of_msgs = prod_of_msgs .* imsg(1:nStates(n),e);
			else
				prod_of_msgs = prod_of_msgs .* imsg(1:nStates(n),e+nEdges);
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
end

%% Decoding
[~,nodeLabels] = max(nodeBel,[],2);

