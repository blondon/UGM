function [nodeBel,edgeBel,logZ,H] = UGM_Infer_ConvexBP(convexity,nodePot,edgePot,edgeStruct,inferFunc,varargin)
%
% Convexified BP inference.
%
% convexity : scales modulus of convexity of entropy (approximation) by
%				inverse scaling of log-potentials
% nodePot : node potentials
% edgePot : edge potentials
% edgeStruct : edge struct
% inferFunc : inference function

assert(nargin >= 5, 'USAGE: UGM_Infer_ConvexBP(convexity,nodePot,edgePot,edgeStruct,inferFunc)')

assert(convexity > 0, sprintf('Convexity must be strictly positive; convexity=%f',convexity));

% Multiply log-potentials by (1/convexity) and rescale
logNodePot = log(nodePot) / convexity;
logEdgePot = log(edgePot) / convexity;
maxLogNodePot = max(logNodePot,[],2);
maxLogEdgePot = max(max(logEdgePot));
for j = 1:size(nodePot,2)
	logNodePot(:,j) = logNodePot(:,j) - maxLogNodePot;
end
for j = 1:size(edgePot, 1)
	for k = 1:size(edgePot,2)
		logEdgePot(j,k,:) = logEdgePot(j,k,:) - maxLogEdgePot;
	end
end
nodePot = exp(logNodePot);
edgePot = exp(logEdgePot);

% Inference
if nargout == 1
	[nodeBel] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});
elseif nargout == 2
	[nodeBel,edgeBel] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});
else
	[nodeBel,edgeBel,logZ,H] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});
end
