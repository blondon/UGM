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

%% Find tree weights
if isscalar(mu) % Weights not provided, construct them using one of the methods below
	% Compute Edge Appearance Probabilities
	if mu == 0
		mu = ones(nEdges,1); % Ordinary BP (not a valid distribution over trees, so not convex)
	elseif mu == 1
		% Generate Random Spanning Trees until all edges are covered
		[nNodes,maxState] = size(nodePot);
		edgeEnds = edgeStruct.edgeEnds;
		
		i = 0;
		edgeAppears = zeros(nEdges,1);
		while 1
			i = i+1;
			edgeAppears = edgeAppears+minSpan(nNodes,[edgeEnds rand(nEdges,1)]);
			if all(edgeAppears > 0)
				break;
			end
		end
		mu = edgeAppears/i;
	elseif mu == 2
		% Compute all spanning trees of the dense graph (not a valid distribution over trees for over graphs)
		mu = ((nNodes-1)/nEdges)*ones(nEdges,1);
	end
end

%% BP
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
