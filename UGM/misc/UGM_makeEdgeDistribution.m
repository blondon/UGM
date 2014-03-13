function mu = UGM_makeEdgeDistribution(edgeStruct,type,maxIters)
%
% Computes the edge distribution for TRBP.
%

if nargin < 2
	type = 1;
end
if nargin < 3
	maxIters = 1000;
end

nNode = edgeStruct.nNodes;
nEdge = edgeStruct.nEdges;

switch(type)
	
	case 0
		% Ordinary BP (not a valid distribution over trees, so not convex)
		mu = ones(nEdge,1);
		
	case 1
		% Generate Random Spanning Trees until all edges are covered
		edgeEnds = edgeStruct.edgeEnds;
		i = 0;
		edgeAppears = zeros(nEdge,1);
		while i < maxIters
			i = i+1
			weights = rand(nEdge,1) .* (1+edgeAppears);
			edgeAppears = edgeAppears+minSpan(nNode,[edgeEnds weights]);
			if all(edgeAppears > 0)
				break;
			end
		end
		% Pretend all uncovered edges got covered
		% (This is not correct; probably should use Laplacian smoothing)
		if any(edgeAppears == 0)
			warning('UGM_makeEdgeDistribution : did not cover %d edges after %d iterations', nnz(edgeAppears==0),maxIters);
			edgeAppears(edgeAppears==0) = 1;
		end
		% Normalize by number of trees
		mu = edgeAppears/i;
		
	case 2
		% Compute all spanning trees of the dense graph (not a valid distribution over trees for over graphs)
		mu = ((nNode-1)/nEdge) * ones(nEdge,1);
		
	case 3
		% Approximate edge distribution of a grid
		mu = zeros(nEdge,1);
		edgeEnds = edgeStruct.edgeEnds;
		for e = 1:nEdge
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);
			% Check for boundary edge
			%TODO
		end
		
end


