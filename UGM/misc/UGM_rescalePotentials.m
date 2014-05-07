function [nodePotSmall,edgePotSmall] = UGM_rescalePotentials(nodePot,edgePot)
%
% Rescales the node/edge potentials by their respective maxima.
%
% Note: This will not affect the beliefs returned by inference, but it will
% affect the partition.

maxNodePot = max(nodePot,[],2);
maxEdgePot = max(max(edgePot));
nodePotSmall = nodePot;
for j = 1:size(nodePot,2)
	nodePotSmall(:,j) = nodePot(:,j)./maxNodePot;
end
edgePotSmall = edgePot;
for j = 1:size(edgePot, 1)
	for k = 1:size(edgePot,2)
		edgePotSmall(j,k,:) = edgePot(j,k,:)./maxEdgePot;
	end
end

