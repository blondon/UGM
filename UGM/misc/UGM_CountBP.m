function [imsg,omsg,convergedStatus] = UGM_CountBP(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,convTol,maximize)

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
imsg_old = imsg;
omsg_old = omsg;

% Precompute aux counting numbers
auxNodeCount = zeros(nNodes,1);
for n = 1:nNodes
	neighbs = UGM_getEdges(n,edgeStruct);
	auxNodeCount(n) = nodeCount(n) + sum(edgeCount(neighbs));
end

% Convergence status
convergedStatus = -1;

% Main loop
for i = 1:edgeStruct.maxIter
	
	% Iterate over nodes, in sequence
	for n = 1:nNodes
		
		% Neighbors of n
		neighbs = UGM_getEdges(n,edgeStruct);
		
		% Incoming
		% log(imsg) = log(
		%	( (\sum_{x_e \ x_n} pot_e(x_e) \prod_{n' : n' \in e\n} omsg_{n',e})^{1/c_e} )^c_e
		% )
		for e = neighbs
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
				tmp = max(bsxfun(@plus, pot_ij, msg') / edgeCount(e), [], 2);
			else
				tmp = logsumexp(bsxfun(@plus, pot_ij, msg')' / edgeCount(e))';
			end
			imsg(:,eidx) = logNormalize(tmp * edgeCount(e));
		end
		
		% Outgoing
		% log(omsg) = log( (\prod_{e' : n \in e'} imsg_{e',n})^{c_e / \hat c_n} / imsg_{e,n} )
		logProd = nodePot(n,:)';
		for e = neighbs
			if n == edgeEnds(e,1)
				logProd = logProd + imsg(:,e);
			else
				logProd = logProd + imsg(:,e+nEdges);
			end
		end
		for e = neighbs
			if n == edgeEnds(e,1)
				omsg(:,e) = logNormalize(logProd .* (edgeCount(e)/auxNodeCount(n)) - imsg(:,e));
			else
				omsg(:,e+nEdges) = logNormalize(logProd .* (edgeCount(e)/auxNodeCount(n)) - imsg(:,e+nEdges));
			end
		end
		
	end
	
	% Check convergence
	if all(abs(imsg(:)-imsg_old(:)) < convTol) && all(abs(omsg(:)-omsg_old(:)) < convTol)
		convergedStatus = i;
		break
	end
	
	% Update old messages
	imsg_old = imsg;
	omsg_old = omsg;
	
end

if convergedStatus == -1
	fprintf('CountBP did not converge after %d iterations\n',edgeStruct.maxIter);
end

% Convert to exponential space
imsg = exp(imsg);
omsg = exp(omsg);

end % END UGM_CountBP



%% "normalize" log-form message
function msg = logNormalize(msg)

msg = msg - max(msg);

end



%% compute log(sum(exp(x)))
function y = logsumexp(x)

maxval = max(x);

y = bsxfun(@plus, log(sum(exp(bsxfun(@minus, x, maxval)))), maxval);

end
