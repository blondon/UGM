function [nodeLabels,mu] = UGM_Decode_TRBP(nodePot,edgePot,edgeStruct,edgeDist)
%
% Approximate decoding based on tree-reweighted beliefe propagation.
% Note: assumes no ties.

if nargin >= 4
	mu = edgeDist;
elseif isfield(edgeStruct,'edgeDist')
	mu = edgeStruct.edgeDist;
else
	mu = 1;
end
[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);

% Find tree weights
if isscalar(mu)
	% Weights not provided, construct them using method specified in mu.
	mu = UGM_makeEdgeDistribution(edgeStruct,mu);
end

%% TRBP
if edgeStruct.useMex
	nodeBel = UGM_Decode_TRBPC(nodePot,edgePot,edgeStruct.edgeEnds,edgeStruct.nStates,edgeStruct.V,edgeStruct.E,int32(edgeStruct.maxIter),mu);
else
	
	maximize = 1;
	new_msg = UGM_TRBP(nodePot,edgePot,edgeStruct,maximize,mu);
	
	[nNodes,maxState] = size(nodePot);
	nEdges = size(edgePot,3);
	edgeEnds = edgeStruct.edgeEnds;
	V = edgeStruct.V;
	E = edgeStruct.E;
	nStates = edgeStruct.nStates;
	
	
	%% Compute nodeBel
	nodeBel = zeros(nNodes,maxState);
	for n = 1:nNodes
        edges = UGM_getEdges(n,edgeStruct);
		prod_of_msgs(1:nStates(n),n) = nodePot(n,1:nStates(n))';
		for e = edges
			if n == edgeEnds(e,2)
				prod_of_msgs(1:nStates(n),n) = prod_of_msgs(1:nStates(n),n) .* (new_msg(1:nStates(n),e).^mu(e));
			else
				prod_of_msgs(1:nStates(n),n) = prod_of_msgs(1:nStates(n),n) .* (new_msg(1:nStates(n),e+nEdges).^mu(e));
			end
		end
		nodeBel(n,1:nStates(n)) = prod_of_msgs(1:nStates(n),n)'./sum(prod_of_msgs(1:nStates(n),n));
	end
end

%% Decoding
[pot nodeLabels] = max(nodeBel,[],2);
