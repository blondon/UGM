function mu = UGM_makeEdgeDistribution(edgeStruct,type,varargin)
%
% Computes the edge distribution for TRBP.
%
% edgeStruct : edge structure
% type : type of distribution
%			0 : P(edge) = 1 uniformly (ordinary BP)
%			1 : (def) Generates random MSTs until all edges covered at least once
%			2 : Computes all spanning trees of the dense graph
%			3 : Approximate edge distribution of a grid (requires [nRows nCols])

if nargin < 2
	type = 1;
end

nNode = edgeStruct.nNodes;
nEdge = edgeStruct.nEdges;

switch(type)
	
	case 0
		% Ordinary BP (not a valid distribution over trees, so not convex)
		mu = ones(nEdge,1);
		
	case 1
		% Generate Random Spanning Trees until all edges are covered
		if isempty(varargin)
			maxIters = 1000;
		else
			maxIters = varargin{1};
		end
		edgeEnds = edgeStruct.edgeEnds;
		i = 0;
		edgeAppears = zeros(nEdge,1);
		while i < maxIters
			i = i+1;
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
		assert(length(varargin)>=1,...
			'USAGE: UGM_makeEdgeDistribution(edgeStruct,3,nRows,nCols) or UGM_makeEdgeDistribution(edgeStruct,3,[nRows nCols])')
		dims = cell2mat(varargin);
		nRows = dims(1);
		nCols = dims(2);
		mu = zeros(nEdge,1);
		edgeEnds = edgeStruct.edgeEnds;
		for e = 1:nEdge
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);
			% Check whether edge is vertical or horizontal
			%   If the edge is vertical, |n1-n2| = 1
			if abs(n1-n2) == 1
				% Vertical edge probability is (nRows+1)/(nRows+nCols)
				mu(e) = (nRows+1) / (nRows+nCols);
			else
				% Horizontal edge probability is (nCols+1)/(nRows+nCols)
				mu(e) = (nCols+1) / (nRows+nCols);
			end
		end
		
	case 4
		% "Uniform 4 Comb" distribution over a grid graph
		assert(length(varargin)>=1,...
			'USAGE: UGM_makeEdgeDistribution(edgeStruct,3,nRows,nCols) or UGM_makeEdgeDistribution(edgeStruct,3,[nRows nCols])')
		dims = cell2mat(varargin);
		nRows = dims(1);
		nCols = dims(2);
		mu = zeros(nEdge,1);
		edgeEnds = edgeStruct.edgeEnds;
		for e = 1:nEdge
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);
			% Check whether edge is vertical or horizontal
			%   If the edge is vertical, |n1-n2| = 1
			if abs(n1-n2) == 1
				% Vertical edge
				if n1 <= nRows || n1 > nRows*(nCols-1)
					% Pr( ver edge on sides ) = 3/4
					mu(e) = 0.75;
				else
					% Pr( ver edge not on sides ) = 1/2
					mu(e) = 0.5;
				end
			else
				% Horizontal edge
				if mod(n1,nRows) == 0 || mod(n1,nRows) == 1
					% Pr( hor edge on sides ) = 3/4
					mu(e) = 0.75;
				else
					% Pr( hor edge not on sides ) = 1/2
					mu(e) = 0.5;
				end
			end
		end

end


