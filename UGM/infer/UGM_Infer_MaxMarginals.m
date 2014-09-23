function [nodeBel] = UGM_Infer_MaxMarginals(nodePot,edgePot,edgeStruct)
% Assumes no ties

if edgeStruct.useMex
	nodeBel = UGM_Decode_LBPC(nodePot,edgePot,edgeStruct.edgeEnds,edgeStruct.nStates,edgeStruct.V,edgeStruct.E,int32(edgeStruct.maxIter));
else
	maximize = 1;
	new_msg = UGM_LoopyBP(nodePot,edgePot,edgeStruct,maximize);
	
	[nNodes,maxState] = size(nodePot);
	nEdges = size(edgePot,3);
	edgeEnds = edgeStruct.edgeEnds;
	nStates = edgeStruct.nStates;

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

% 	if nargout > 1
% 		% Compute edge beliefs
% 		edgeBel = zeros(maxState,maxState,nEdges);
% 		for e = 1:nEdges
% 			n1 = edgeEnds(e,1);
% 			n2 = edgeEnds(e,2);
% 			belN1 = nodeBel(n1,1:nStates(n1))'./new_msg(1:nStates(n1),e+nEdges);
% 			belN2 = nodeBel(n2,1:nStates(n2))'./new_msg(1:nStates(n2),e);
% 			b1=repmat(belN1,1,nStates(n2));
% 			b2=repmat(belN2',nStates(n1),1);
% 			eb = b1.*b2.*edgePot(1:nStates(n1),1:nStates(n2),e);
% 			edgeBel(1:nStates(n1),1:nStates(n2),e) = eb./sum(eb(:));
% 		end
% 	end
end

