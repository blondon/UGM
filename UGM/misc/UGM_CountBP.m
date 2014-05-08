function [msg_i,msg_o] = UGM_CountBP(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,maximize)

[nNodes,maxState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
nStates = double(edgeStruct.nStates);

% Initialize
msg_i = zeros(maxState,nEdges*2);
msg_o = zeros(maxState,nEdges*2);
for e = 1:nEdges
	n1 = edgeEnds(e,1);
	n2 = edgeEnds(e,2);
	msg_i(1:nStates(n1),e) = 1/nStates(n1); % Message e -> n1
	msg_i(1:nStates(n2),e+nEdges) = 1/nStates(n2); % Message e -> n2
	msg_o(1:nStates(n1),e) = 1/nStates(n1); % Message n1 -> e
	msg_o(1:nStates(n2),e+nEdges) = 1/nStates(n2); % Message n2 -> e
	edgePot(1:nStates(n1),1:nStates(n2),e) = edgePot(1:nStates(n1),1:nStates(n2),e).^(1/edgeCount(e));
end
old_msg_i = msg_i;
old_msg_o = msg_o;

% Main loop
for i = 1:edgeStruct.maxIter
	% Temp messages
	tmp_i = zeros(maxState,nNodes);
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
		
		% exponent variables
		q1 = (1 - nodeCount(n1)) / length(UGM_getEdges(n1,edgeStruct));
		q2 = (1 - nodeCount(n2)) / length(UGM_getEdges(n2,edgeStruct));
		d1 = edgeCount(e) - q1 + 1;
		d2 = edgeCount(e) - q2 + 1;
		
		% Incoming
		%  n1
		newm = tmp_i(1:nStates(n1),e).^(edgeCount(e)/d1) .* ...
			   tmp_o(1:nStates(n1),e).^((q1-edgeCount(e))/d1);
		newm(isnan(newm)|isinf(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
		end
		msg_i(1:nStates(n1),e) = newm;
		%  n2
		newm = tmp_i(1:nStates(n2),e+nEdges).^(edgeCount(e)/d2) .* ...
			   tmp_o(1:nStates(n2),e+nEdges).^((q2-edgeCount(e))/d2);
		newm(isnan(newm)|isinf(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
		end
		msg_i(1:nStates(n2),e+nEdges) = newm;

		% Outgoing
		%  n1
		newm = tmp_i(1:nStates(n1),e).^((q1-1)/d1) .* ...
			   tmp_o(1:nStates(n1),e).^(1/d1);
		newm(isnan(newm)|isinf(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
		end
		msg_o(1:nStates(n1),e) = newm;
		%  n2
		newm = tmp_i(1:nStates(n2),e+nEdges).^((q2-1)/d2) .* ...
			   tmp_o(1:nStates(n2),e+nEdges).^(1/d2);
		newm(isnan(newm)|isinf(newm)) = 0;
		z = sum(newm);
		if z > 0
			newm = newm ./ z;
		end
		msg_o(1:nStates(n2),e+nEdges) = newm;
	end
	
% 	% Check for NaNs
% 	[min(msg_i(:)) max(msg_i(:)) min(msg_o(:)) max(msg_o(:))]
% 	if any(isnan(msg_i(:))) || any(isnan(msg_o(:)))
% 		error('Found NaN values\n')
% 	end

	% Damping
	if momentum < 1
		msg_i = old_msg_i.^(1-momentum) .* msg_i.^momentum;
		msg_o = old_msg_o.^(1-momentum) .* msg_o.^momentum;
	end
	
	% Check convergence
	if all(abs(msg_i(:)-old_msg_i(:)) < convTol) ...
	   && all(abs(msg_o(:)-old_msg_o(:)) < convTol)
		break
	end
	
	% Swap old/new outgoing messages
	old_msg_i = msg_i;
	old_msg_o = msg_o;
end
if i == edgeStruct.maxIter
	fprintf('CountBP did not converge after %d iterations\n',edgeStruct.maxIter);
end
% fprintf('CountBP stopped after %d iterations\n',i);
