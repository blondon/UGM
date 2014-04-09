function [nodeCount,edgeCount] = UGM_TRBPCounts(edgeStruct,edgeDist)
%
% Computes the counting numbers for TRBP.
%
% edgeStruct : edge structure
% edgeDist : edge distribution; i.e., P(e)
%				(optional; if not supplied, must be field in edgeStruct)
%
% nodeCount : (nNodes x 1) vector of node counting numbers, where
%				nodeCount(n) = 1 - sum_{e : n \in e} P(e)
% edgeCount : (nEdges x 1) vector of edge counting numbers, where
%				nodeCount(e) = P(e)

assert(nargin == 2 || isfield(edgeStruct,'edgeDist'),...
	'USAGE: UGM_TRBPCounts(edgeStruct), with edgeStruct.edgeDist, or UGM_TRBPCounts(edgeStruct,edgeDist).');

if nargin < 2
	edgeDist = edgeStruct.edgeDist;
end

nodeCount = ones(edgeStruct.nNodes,1);
edgeCount = edgeDist;
for e = 1:edgeStruct.nEdges
	n1 = edgeStruct.edgeEnds(e,1);
	n2 = edgeStruct.edgeEnds(e,2);
	nodeCount(n1) = nodeCount(n1) - edgeDist(e);
	nodeCount(n2) = nodeCount(n2) - edgeDist(e);
end
