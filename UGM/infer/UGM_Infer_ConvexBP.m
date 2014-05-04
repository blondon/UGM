function [nodeBel,edgeBel,logZ] = UGM_Infer_ConvexBP(convexity,nodePot,edgePot,edgeStruct,inferFunc,varargin)
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

if nargout == 1
	[nodeBel] = inferFunc(nodePot.^(1/convexity),edgePot.^(1/convexity),edgeStruct,varargin{:});
elseif nargout == 2
	[nodeBel,edgeBel] = inferFunc(nodePot.^(1/convexity),edgePot.^(1/convexity),edgeStruct,varargin{:});
else
	[nodeBel,edgeBel,logZ] = inferFunc(nodePot.^(1/convexity),edgePot.^(1/convexity),edgeStruct,varargin{:});
end
