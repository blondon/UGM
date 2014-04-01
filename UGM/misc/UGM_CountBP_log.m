function [msg_i,msg_o] = UGM_CountBP_log(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,maximize)

[nNodes,maxState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
nStates = double(edgeStruct.nStates);

% Convergence tolerance
if isfield(edgeStruct,'convTol')
	convTol = edgeStruct.convTol;
else
	convTol = 1e-4;
end

% Initialize
lnodePot = log(nodePot);
lmsg_i = zeros(maxState,nEdges*2);
lmsg_o = zeros(maxState,nEdges*2);
lold_msg_i = zeros(maxState,nEdges*2);
lold_msg_o = zeros(maxState,nEdges*2);
for e = 1:nEdges
	n1 = edgeEnds(e,1);
	n2 = edgeEnds(e,2);
	lmsg_i(1:nStates(n1),e) = -log(nStates(n1)); % Message e -> n1
	lmsg_i(1:nStates(n2),e+nEdges) = -log(nStates(n2)); % Message e -> n2
	lmsg_o(1:nStates(n1),e) = -log(nStates(n1)); % Message n1 -> e
	lmsg_o(1:nStates(n2),e+nEdges) = -log(nStates(n2)); % Message n2 -> e
	edgePot(1:nStates(n1),1:nStates(n2),e) = edgePot(1:nStates(n1),1:nStates(n2),e).^(1/edgeCount(e));
end

% Main loop
for i = 1:edgeStruct.maxIter
	% Temp messages
	ltmp_i = zeros(maxState,nNodes);
	ltmp_o = zeros(maxState,nEdges*2);
	for e = 1:nEdges
		n1 = edgeEnds(e,1);
		n2 = edgeEnds(e,2);
		
		% Incoming
		pot = edgePot(1:nStates(n1),1:nStates(n2),e);
		z = min(lmsg_o(1:nStates(n2),e+nEdges));
		ltmp_i(1:nStates(n1),e) = log(pot * exp(lmsg_o(1:nStates(n2),e+nEdges)-z)) + z;
		z = min(lmsg_o(1:nStates(n1),e));
		ltmp_i(1:nStates(n2),e+nEdges) = log(pot' * exp(lmsg_o(1:nStates(n1),e)-z)) + z;
		
		% Outgoing
		%  n1
		tmp = lnodePot(n1,1:nStates(n1))';
		for e1 = UGM_getEdges(n1,edgeStruct)
			if e1 ~= e
				if n1 == edgeEnds(e1,1)
					tmp = tmp + lmsg_i(1:nStates(n1),e1);
				else
					tmp = tmp + lmsg_i(1:nStates(n1),e1+nEdges);
				end
			end
		end
		ltmp_o(1:nStates(n1),e) = tmp;
		%  n2
		tmp = lnodePot(n2,1:nStates(n2))';
		for e2 = UGM_getEdges(n2,edgeStruct)
			if e2 ~= e
				if n2 == edgeEnds(e2,1)
					tmp = tmp + lmsg_i(1:nStates(n2),e2);
				else
					tmp = tmp + lmsg_i(1:nStates(n2),e2+nEdges);
				end
			end
		end
		ltmp_o(1:nStates(n2),e+nEdges) = tmp;
	end
	
	% New messages
	for e = 1:nEdges
		n1 = edgeEnds(e,1);
		n2 = edgeEnds(e,2);
		
		% exponent variables
		q1 = (1 - nodeCount(n1)) / length(UGM_getEdges(n1,edgeStruct));
		q2 = (1 - nodeCount(n2)) / length(UGM_getEdges(n2,edgeStruct));
		d1 = edgeCount(e) - q1 + 1;
		d2 = edgeCount(e) - q2 + 1;
		
		% Incoming
		%  n1
		newm = ltmp_i(1:nStates(n1),e).*(edgeCount(e)/d1) + ...
			   ltmp_o(1:nStates(n1),e).*((q1-edgeCount(e))/d1);
		lmsg_i(1:nStates(n1),e) = newm - log(sum(exp(newm)));
		%  n2
		newm = ltmp_i(1:nStates(n2),e+nEdges).*(edgeCount(e)/d2) + ...
			   ltmp_o(1:nStates(n2),e+nEdges).*((q2-edgeCount(e))/d2);
		lmsg_i(1:nStates(n2),e+nEdges) = newm - log(sum(exp(newm)));

		% Outgoing
		newm = ltmp_i(1:nStates(n1),e).*((q1-1)/d1) + ...
			   ltmp_o(1:nStates(n1),e).*(1/d1);
		lmsg_o(1:nStates(n1),e) = newm - log(sum(exp(newm)));
		newm = ltmp_i(1:nStates(n2),e+nEdges).*((q2-1)/d2) + ...
			   ltmp_o(1:nStates(n2),e+nEdges).*(1/d2);
		lmsg_o(1:nStates(n2),e+nEdges) = newm - log(sum(exp(newm)));
	end
% 	
% 	% Check for NaNs
% 	if any(isnan(lmsg_i(:))) || any(isnan(lmsg_o(:)))
% 		break
% 	end
	
	% Check convergence
	if all(abs(exp(lmsg_i(:))-exp(lold_msg_i(:))) < convTol) ...
	   && all(abs(exp(lmsg_o(:))-exp(lold_msg_o(:))) < convTol)
		break
	end
	
	% Swap old/new outgoing messages
	lold_msg_i = lmsg_i;
	lold_msg_o = lmsg_o;
end
if i == edgeStruct.maxIter
	fprintf('CountBP did not converge after %d iterations\n',edgeStruct.maxIter);
end
fprintf('CountBP stopped after %d iterations\n',i);

msg_i = exp(lmsg_i);
msg_o = exp(lmsg_o);
% for e = 1:nEdges
% 	n1 = edgeEnds(e,1);
% 	n2 = edgeEnds(e,2);
% 	msg_i(1:nStates(n1),e) = msg_i(1:nStates(n1),e) ./ sum(msg_i(1:nStates(n1),e));
% 	msg_i(1:nStates(n2),e+nEdges) = msg_i(1:nStates(n2),e+nEdges) ./ sum(msg_i(1:nStates(n2),e+nEdges));
% 	msg_o(1:nStates(n1),e) = msg_o(1:nStates(n1),e) ./ sum(msg_o(1:nStates(n1),e));
% 	msg_o(1:nStates(n2),e+nEdges) = msg_o(1:nStates(n2),e+nEdges) ./ sum(msg_o(1:nStates(n2),e+nEdges));
% end
