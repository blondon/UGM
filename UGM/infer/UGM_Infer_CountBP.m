function [nodeBel,edgeBel,logZ] = UGM_Infer_CountBP(nodePot,edgePot,edgeStruct,nodeCount,edgeCount)

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
    [nodeBel,edgeBel,logZ] = UGM_Infer_CountBPC(...
		nodePot,edgePot,nodeCount,edgeCount,...
		edgeStruct.edgeEnds,edgeStruct.nStates,edgeStruct.V,edgeStruct.E,...
		int32(edgeStruct.maxIter),momentum,convTol);
else
    [nodeBel, edgeBel, logZ] = Infer_CountBP(nodePot,edgePot,edgeStruct,...
		nodeCount,edgeCount,momentum,convTol);
end

end

%% Non-mex version
function [nodeBel,edgeBel,logZ] = Infer_CountBP(nodePot,edgePot,edgeStruct,nodeCount,edgeCount,momentum,convTol)

nNodes = edgeStruct.nNodes;
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;
maxStates = max(nStates);
V = edgeStruct.V;
E = edgeStruct.E;

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
    nodeBel(n,1:nStates(n)) = prod_of_msgs' ./ sum(prod_of_msgs);
end

if nargout > 1
    % Compute edge beliefs
    edgeBel = zeros(maxStates,maxStates,nEdges);
	for e = 1:nEdges
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);
		eb = edgePot(1:nStates(n1),1:nStates(n2),e) .* ...
			 (msg_o(1:nStates(n1),e) * msg_o(1:nStates(n2),e+nEdges)');
		edgeBel(1:nStates(n1),1:nStates(n2),e) = eb ./ sum(eb(:));
	end
end

if nargout > 2
    % Compute Bethe free energy
    Energy1 = 0; Energy2 = 0; Entropy1 = 0; Entropy2 = 0;
    nodeBel = nodeBel+eps;
    edgeBel = edgeBel+eps;
    for n = 1:nNodes
        edges = E(V(n):V(n+1)-1);
        nNbrs = length(edges);

        % Node Entropy (can get divide by zero if beliefs at 0)
        Entropy1 = Entropy1 + (nNbrs-1)*sum(nodeBel(n,1:nStates(n)).*log(nodeBel(n,1:nStates(n))));

        % Node Energy
        Energy1 = Energy1 - sum(nodeBel(n,1:nStates(n)).*log(nodePot(n,1:nStates(n))));
    end
    for e = 1:nEdges
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);

        % Pairwise Entropy (can get divide by zero if beliefs at 0)
        eb = edgeBel(1:nStates(n1),1:nStates(n2),e);
        Entropy2 = Entropy2 - sum(eb(:).*log(eb(:)));

        % Pairwise Energy
        ep = edgePot(1:nStates(n1),1:nStates(n2),e);
        Energy2 = Energy2 - sum(eb(:).*log(ep(:)));
    end
    F = (Energy1+Energy2) - (Entropy1+Entropy2);
    logZ = -F;
end

end