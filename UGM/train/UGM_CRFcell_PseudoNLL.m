function [nll,g] = UGM_CRFcell_PseudoNLL(w,examples,varargin)
% [nll,g] = UGM_CRFcell_PseudoNLL(w,examples,varargin)
%
%   calls UGM_CRF_PseudoNLL for each set of examples in the cell array examples,
%   and sums the result

nll = 0;
g = zeros(size(w));
for i = 1:length(examples)
   [nll_sub,g_sub] = UGM_CRF_PseudoNLL(w,examples{i}.Xnode,examples{i}.Xedge,examples{i}.Y',examples{i}.nodeMap,examples{i}.edgeMap,examples{i}.edgeStruct,varargin{:});
    nll = nll + nll_sub;
    g = g + g_sub;
end