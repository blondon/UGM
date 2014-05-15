function [msg] = UGM_CountBP_new(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,maximize)

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
	edgePot(:,:,e) = edgePot(:,:,e).^(1/edgeCount(e));
end

% Init messages
msg = (1/nState) * ones(nState,nEdges*2); % n -> e
msg_old = msg;
itmp = zeros(nState,nEdges*2);
otmp = zeros(nState,nEdges*2);

% Main loop
for i = 1:edgeStruct.maxIter

	% Iterate over nodes, in sequence. (Crucial to convergence!!!)
	for n = 1:nNodes
		
		% Neighbors of n
		neighbs = UGM_getEdges(n,edgeStruct);
		
		% Update temp variables
		for e = neighbs
			% Incoming
			if n == edgeEnds(e,1)
				eidx = e;
				pot_ij = edgePot(:,:,e);
				m = msg(:,e+nEdges);
			else
				eidx = e + nEdges;
				pot_ij = edgePot(:,:,e)';
				m = msg(:,e);
			end
			if maximize
				itmp(:,eidx) = max(pot_ij .* repmat(m',nState,1), [], 2);
			else
				itmp(:,eidx) = pot_ij * m;
			end
			% Outgoing
			prod = nodePot(n,:)';
			for e_ = neighbs
				if e_ ~= e
					if n == edgeEnds(e_,1)
						prod = prod .* msg(:,e_+nEdges);
					else
						prod = prod .* msg(:,e_);
					end
				end
			end
			otmp(:,eidx) = prod;
		end
		
		% Update messages
		for e = neighbs
			if n == edgeEnds(e,1)
				eidx1 = e;
				eidx2 = e + nEdges;
			else
				eidx1 = e + nEdges;
				eidx2 = e;
			end
			% Incoming
			msg(:,eidx2) = normalizeMessage(itmp(:,eidx2).^iexp(eidx2) .* otmp(:,eidx2).^(oexp(eidx2)-1));
			% Outgoing
			msg(:,eidx1) = normalizeMessage(itmp(:,eidx1).^(iexp(eidx1)-1) .* otmp(:,eidx1).^oexp(eidx1));
		end
		
	end

	% Damping
	if momentum < 1
		msg = msg_old.^(1-momentum) .* msg.^momentum;
	end
	
	% Check convergence
	fprintf('sum-abs-diff = %f, max-abs-diff = %f\n', ...
		sum(abs(msg(:)-msg_old(:))), max(abs(msg(:)-msg_old(:))));
	if all(abs(msg(:)-msg_old(:)) < convTol)
		break
	end
	
	% Update old messages
	msg_old = msg;
	
end

if i == edgeStruct.maxIter
	fprintf('CountBP did not converge after %d iterations\n',edgeStruct.maxIter);
end
fprintf('CountBP stopped after %d iterations\n',i);

end % END UGM_CountBP

%% Normalize a message while removing Inf and NaN
function m = normalizeMessage(m)

m(~isfinite(m)) = 0;
z = sum(m);
if z > 0
	m = m ./ z;
end

end % END normalizeMessage
