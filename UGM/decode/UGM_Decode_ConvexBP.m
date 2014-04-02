function  [nodeLabels] = UGM_Decode_ConvexBP(convexity,nodePot,edgePot,edgeStruct,inferFunc,varargin)
%
% Decode convexified BP inference.
%
% convexity : scales modulus of convexity of entropy (approximation) by
%				inverse scaling of log-potentials
% nodePot : node potentials
% edgePot : edge potentials
% edgeStruct : edge struct
% inferFunc : inference function

assert(nargin >= 5, 'USAGE: UGM_Decode_ConvexBP(convexity,nodePot,edgePot,edgeStruct,inferFunc)')

assert(convexity > 0, 'In UGM_Decode_ConvexBP: convexity must be strictly positive.')

nodeBel = UGM_Infer_ConvexBP(convexity,nodePot,edgePot,edgeStruct,inferFunc,varargin{:});
[~,nodeLabels] = max(nodeBel,[],2);
