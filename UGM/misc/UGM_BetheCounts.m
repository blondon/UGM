function [nodeCount,edgeCount] = UGM_BetheCounts(edgeStruct)
%
% Computes the counting numbers for the Bethe approximation.
%
% edgeStruct : edge structure
%
% nodeCount : (nNodes x 1) vector of node counting numbers, where
%				nodeCount(n) = 1 - deg(n)
% edgeCount : (nEdges x 1) vector of edge counting numbers, where
%				nodeCount(e) = 1

nodeCount = ones(edgeStruct.nNodes,1);
edgeCount = ones(edgeStruct.nEdges,1);

for e = 1:edgeStruct.nEdges
	n1 = edgeStruct.edgeEnds(e,1);
	n2 = edgeStruct.edgeEnds(e,2);
	nodeCount(n1) = nodeCount(n1) - 1;
	nodeCount(n2) = nodeCount(n2) - 1;
end

