function [imsg,omsg] = UGM_CountBP_schwing(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,maximize)

[nNodes,nState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;

% For simplicity, all variables must have same number of states
assert(all(edgeStruct.nStates == nState), 'UGM_CountBP: All variables must have same number of states.')

% Convert potentials to log space
nodePot = log(nodePot);
edgePot = log(edgePot);

% Init messages
imsg = zeros(nState,nEdges*2); % incoming: e -> n
omsg = zeros(nState,nEdges*2); % outgoing: n -> e
itmp = zeros(nState,nEdges*2); % incoming: e -> n
otmp = zeros(nState,nEdges*2); % outgoing: n -> e
imsg_old = imsg;
omsg_old = omsg;

% Precompute aux counting numbers
auxNodeCount = zeros(nNodes,1);
for n = 1:nNodes
	neighbs = UGM_getEdges(n,edgeStruct);
	auxNodeCount(n) = nodeCount(n) + sum(edgeCount(neighbs));
end

% Main loop
for i = 1:edgeStruct.maxIter
	
	% Iterate over nodes, in sequence
	for n = 1:nNodes
		
		% Neighbors of n
		neighbs = UGM_getEdges(n,edgeStruct);
		
		% Update temp variables
		for e = neighbs
			% Incoming
			% log(imsg) = log(
			%	( (\sum_{x_e \ x_n} pot_e(x_e) \prod_{n' : n' \in e\n} omsg_{n',e})^{1/c_e} )^c_e
			% )
			if n == edgeEnds(e,1)
				eidx = e;
				pot_ij = edgePot(:,:,e);
				msg = omsg(:,e+nEdges);
			else
				eidx = e + nEdges;
				pot_ij = edgePot(:,:,e)';
				msg = omsg(:,e);
			end
			if maximize
				itmp(:,eidx) = max((pot_ij + repmat(msg',nState,1)) / edgeCount(e), [], 2);
			else
				itmp(:,eidx) = logsumexp(bsxfun(@plus, pot_ij, msg')' / edgeCount(e))';
			end
			itmp(:,eidx) = itmp(:,eidx) * edgeCount(e);
			% Outgoing
			% log(omsg) = log( (\prod_{e' : n \in e'} imsg_{e',n})^{c_e / \hat c_n} / imsg_{e,n} )
			logProd = nodePot(n,:)';
			for e_ = neighbs
				if n == edgeEnds(e_,1)
					logProd = logProd + imsg(:,e_);
				else
					logProd = logProd + imsg(:,e_+nEdges);
				end
			end
			otmp(:,eidx) = logProd * (edgeCount(e)/auxNodeCount(n)) - itmp(:,eidx);
		end
		
		% Update messages
		for e = neighbs
			if n == edgeEnds(e,1)
				eidx = e;
			else
				eidx = e + nEdges;
			end
			imsg(:,eidx) = logNormalizeMessage(itmp(:,eidx));
			omsg(:,eidx) = logNormalizeMessage(otmp(:,eidx));
		end
		
	end
	
	% Damping
	if momentum < 1
		imsg = imsg_old.*(1-momentum) + imsg.*momentum;
		omsg = omsg_old.*(1-momentum) + omsg.*momentum;
	end
	
	% Check convergence
	fprintf('sum-abs-diff = %f, max-abs-diff = %f\n', ...
		sum([abs(imsg(:)-imsg_old(:));abs(omsg(:)-omsg_old(:))]), ...
		max([abs(imsg(:)-imsg_old(:));abs(omsg(:)-omsg_old(:))]));
	if all(abs(imsg(:)-imsg_old(:)) < convTol) && all(abs(omsg(:)-omsg_old(:)) < convTol)
		break
	end
	
	% Update old messages
	imsg_old = imsg;
	omsg_old = omsg;
	
end

if i == edgeStruct.maxIter
	fprintf('CountBP did not converge after %d iterations\n',edgeStruct.maxIter);
end
fprintf('CountBP stopped after %d iterations\n',i);

% Convert to exponential space
imsg = exp(imsg);
omsg = exp(omsg);

end % END UGM_CountBP



%% "normalize" log-form message
function msg = logNormalizeMessage(msg)

msg = msg - max(msg);

end



%% compute log(sum(exp(x)))
function y = logsumexp(x)

maxval = max(x);

y = bsxfun(@plus, log(sum(exp(bsxfun(@minus, x, maxval)))), maxval);

end
