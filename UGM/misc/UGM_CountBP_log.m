function [imsg,omsg] = UGM_CountBP_log(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,maximize)

[nNodes,nState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;

% For simplicity, all variables must have same number of states.
assert(all(edgeStruct.nStates == nState), 'UGM_CountBP: All variables must have same number of states.')

% Precompute exponents/potentials
iexp = zeros(nEdges*2,1);
oexp = zeros(nEdges*2,1);
for e = 1:nEdges
	n1 = edgeEnds(e,1);
	n2 = edgeEnds(e,2);
	q = (1 - nodeCount(n1)) / length(UGM_getEdges(n1,edgeStruct));
	iexp(e) = edgeCount(e) / (edgeCount(e) - q + 1);
	oexp(e) = 1 / (edgeCount(e) - q + 1);
	q = (1 - nodeCount(n2)) / length(UGM_getEdges(n2,edgeStruct));
	iexp(e+nEdges) = edgeCount(e) / (edgeCount(e) - q + 1);
	oexp(e+nEdges) = 1 / (edgeCount(e) - q + 1);
	% Raise edge potentials to power (1/edgeCount)
	edgePot(:,:,e) = log(edgePot(:,:,e)).*(1/edgeCount(e));
end

% Init messages
imsg = (1/nState) * ones(nState,nEdges*2); % incoming: e -> n
omsg = (1/nState) * ones(nState,nEdges*2); % outgoing: n -> e
imsg_old = imsg;
omsg_old = omsg;
itmp = zeros(nState,nEdges*2);
otmp = zeros(nState,nEdges*2);

% Main loop
for i = 1:edgeStruct.maxIter

	% Iterate over nodes, in sequence
	for n = 1:nNodes
		
		% Neighbors of n
		neighbs = UGM_getEdges(n,edgeStruct);
		
		% Update temp variables
		for e = neighbs
			% Incoming
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
				itmp(:,eidx) = max(pot_ij + repmat(msg',nState,1), [], 2);
			else
				itmp(:,eidx) = logsumexp(bsxfun(@plus, pot_ij, msg')')';
                % should be equivalent to log(exp(pot_ij)*exp(msg))
			end
			% Outgoing
			prod = nodePot(n,:)';
			for e_ = neighbs
				if e_ ~= e
					if n == edgeEnds(e_,1)
						prod = prod + imsg(:,e_);
					else
						prod = prod + imsg(:,e_+nEdges);
					end
				end
			end
			otmp(:,eidx) = prod;
		end
		
		% Update messages
		for e = neighbs
			if n == edgeEnds(e,1)
				eidx = e;
			else
				eidx = e + nEdges;
			end
			imsg(:,eidx) = logNormalizeMessage(itmp(:,eidx).*iexp(eidx) + otmp(:,eidx).*(oexp(eidx)-1));
			omsg(:,eidx) = logNormalizeMessage(itmp(:,eidx).*(iexp(eidx)-1) + otmp(:,eidx).*oexp(eidx));
		end
		
	end

	% Damping
	if momentum < 1
		imsg = imsg_old.*(1-momentum) + imsg.*momentum;
		omsg = omsg_old.*(1-momentum) + omsg.*momentum;
    end
    
    plot([imsg(:); omsg(:)]);
    drawnow;
	
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

end % END UGM_CountBP


%% Normalize a message while removing Inf and NaN
function msg = normalizeMessage(msg)

msg(~isfinite(msg)) = 0;
z = sum(msg);
if z > 0
	msg = msg ./ z;
end

end % END normalizeMessage



%% "normalize" log-form message
function msg = logNormalizeMessage(msg)

msg = msg - max(msg);

end



%% compute log(sum(exp(x)))
function y = logsumexp(x)

maxval = max(x);

y = bsxfun(@plus, log(sum(exp(bsxfun(@minus, x, maxval)))), maxval);

end
