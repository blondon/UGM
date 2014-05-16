function [msg_i,msg_o] = UGM_CountBP(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,maximize)

[nNodes,maxState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
nStates = double(edgeStruct.nStates);

% Initialize
msg_i = zeros(maxState,nEdges*2);
msg_o = zeros(maxState,nEdges*2);
exp_i = zeros(nEdges*2,1);
exp_o = zeros(nEdges*2,1);
for e = 1:nEdges
	n1 = edgeEnds(e,1);
	n2 = edgeEnds(e,2);
	% Init messages
	msg_i(1:nStates(n1),e) = 1/nStates(n1); % Message e -> n1
	msg_i(1:nStates(n2),e+nEdges) = 1/nStates(n2); % Message e -> n2
	msg_o(1:nStates(n1),e) = 1/nStates(n1); % Message n1 -> e
	msg_o(1:nStates(n2),e+nEdges) = 1/nStates(n2); % Message n2 -> e
	edgePot(1:nStates(n1),1:nStates(n2),e) = edgePot(1:nStates(n1),1:nStates(n2),e).^(1/edgeCount(e));
	% Init exponents
	q = (1 - nodeCount(n1)) / length(UGM_getEdges(n1,edgeStruct));
	exp_i(e) = edgeCount(e) / (edgeCount(e) - q + 1);
	exp_o(e) = 1 / (edgeCount(e) - q + 1);
	q = (1 - nodeCount(n2)) / length(UGM_getEdges(n2,edgeStruct));
	exp_i(e+nEdges) = edgeCount(e) / (edgeCount(e) - q + 1);
	exp_o(e+nEdges) = 1 / (edgeCount(e) - q + 1);
end
% Init old messages
old_msg_i = msg_i;
old_msg_o = msg_o;

% Main loop
for i = 1:edgeStruct.maxIter
	% Temp messages
	tmp_i = zeros(maxState,nEdges*2);
	tmp_o = zeros(maxState,nEdges*2);
	for e = 1:nEdges
		n1 = edgeEnds(e,1);
		n2 = edgeEnds(e,2);
		
		% Incoming
		pot = edgePot(1:nStates(n1),1:nStates(n2),e);
		if maximize
			% Max-product
			tmp_i(1:nStates(n1),e) = max(pot.*repmat(msg_o(1:nStates(n2),e+nEdges)',nStates(n1),1), [], 2);
			tmp_i(1:nStates(n2),e+nEdges) = max(pot'.*repmat(msg_o(1:nStates(n1),e)',nStates(n2),1), [], 2);
		else
			% Sum-product
			tmp_i(1:nStates(n1),e) = pot * msg_o(1:nStates(n2),e+nEdges);
			tmp_i(1:nStates(n2),e+nEdges) = pot' * msg_o(1:nStates(n1),e);
		end
		
		% Outgoing
		%  n1
		tmp = nodePot(n1,1:nStates(n1))';
		for e1 = UGM_getEdges(n1,edgeStruct)
			if e1 ~= e
				if n1 == edgeEnds(e1,1)
					tmp = tmp .* msg_i(1:nStates(n1),e1);
				else
					tmp = tmp .* msg_i(1:nStates(n1),e1+nEdges);
				end
			end
		end
		tmp_o(1:nStates(n1),e) = tmp;
		%  n2
		tmp = nodePot(n2,1:nStates(n2))';
		for e2 = UGM_getEdges(n2,edgeStruct)
			if e2 ~= e
				if n2 == edgeEnds(e2,1)
					tmp = tmp .* msg_i(1:nStates(n2),e2);
				else
					tmp = tmp .* msg_i(1:nStates(n2),e2+nEdges);
				end
			end
		end
		tmp_o(1:nStates(n2),e+nEdges) = tmp;
	end
	
	% New messages
	for e = 1:nEdges
		n1 = edgeEnds(e,1);
		n2 = edgeEnds(e,2);
		
		% Incoming
		%  n1
		newm = tmp_i(1:nStates(n1),e).^exp_i(e) .* ...
			   tmp_o(1:nStates(n1),e).^(exp_o(e)-1);
		newm(~isfinite(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
% 		else
% 			newm = 1 / nStates(n1);
		end
		msg_i(1:nStates(n1),e) = newm;
		%  n2
		newm = tmp_i(1:nStates(n2),e+nEdges).^exp_i(e+nEdges) .* ...
			   tmp_o(1:nStates(n2),e+nEdges).^(exp_o(e+nEdges)-1);
		newm(~isfinite(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
% 		else
% 			newm = 1 / nStates(n2);
		end
		msg_i(1:nStates(n2),e+nEdges) = newm;

		% Outgoing
		%  n1
		newm = tmp_i(1:nStates(n1),e).^(exp_i(e)-1) .* ...
			   tmp_o(1:nStates(n1),e).^exp_o(e);
		newm(~isfinite(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
% 		else
% 			newm = 1 / nStates(n1);
		end
		msg_o(1:nStates(n1),e) = newm;
		%  n2
		newm = tmp_i(1:nStates(n2),e+nEdges).^(exp_i(e+nEdges)-1) .* ...
			   tmp_o(1:nStates(n2),e+nEdges).^exp_o(e+nEdges);
		newm(~isfinite(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
% 		else
% 			newm = 1 / nStates(n2);
		end
		msg_o(1:nStates(n2),e+nEdges) = newm;
	end

	% Damping
	if momentum < 1
		msg_i = old_msg_i.^(1-momentum) .* msg_i.^momentum;
		msg_o = old_msg_o.^(1-momentum) .* msg_o.^momentum;
	end
	
	% Check convergence
	fprintf('sum-abs-diff = %f, max-abs-diff = %f\n', ...
		sum([abs(msg_i(:)-old_msg_i(:));abs(msg_o(:)-old_msg_o(:))]), ...
		max([abs(msg_i(:)-old_msg_i(:));abs(msg_o(:)-old_msg_o(:))]));
	if all(abs(msg_i(:)-old_msg_i(:)) < convTol) && all(abs(msg_o(:)-old_msg_o(:)) < convTol)
		break
	end
% 	% Sum-abs-diff criteria not consistent with mex version
% 	if (abs(msg_i(:)-old_msg_i(:)) + abs(msg_o(:)-old_msg_o(:))) < convTol
% 		break
% 	end
	
	% Update old messages
	old_msg_i = msg_i;
	old_msg_o = msg_o;
end

if i == edgeStruct.maxIter
	fprintf('CountBP did not converge after %d iterations\n',edgeStruct.maxIter);
end
fprintf('CountBP stopped after %d iterations\n',i);
