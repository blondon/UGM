function [bel,imsg,omsg] = UGM_CountBP(nodePot,edgePot,nodeCount,edgeCount,auxCount,edgeStruct,momentum,convTol,maximize)

[nNodes,nStates] = size(nodePot);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;

% For now, only supports variables having same number of states
assert(all(edgeStruct.nStates==nStates), 'UGM_CountBP : all variables must have same number of states');

% Convert counting numbers
c_a = edgeCount;
for e = 1:nEdges
	c_a = c_a - auxCount(e);
	c_a = c_a - auxCount(e+nEdges);
end
c_i = zeros(nNodes,1);
for n = 1:nNodes
	c_i(n) = nodeCount(n);
	for e = UGM_getEdges(n,edgeStruct)
		if n == edgeEnds(e,1)
			c_i(n) = c_i(n) + auxCount(e) + c_a(e);
		else
			c_i(n) = c_i(n) + auxCount(e+nEdges) + c_a(e);
		end
	end
end
c_ia = [c_a ; c_a] + auxCount;

% Initialize
bel = zeros(nStates,nNodes); % node beliefs
bel_old = bel;
imsg = zeros(nStates,nEdges*2); % Incoming message e -> n
omsg = ones(nStates,nStates,nEdges*2); % Outgoing message n1 -> e, (state1,state2,e)

% Main loop
for t = 1:edgeStruct.maxIter
	
	for n = 1:nNodes
		
		% Incoming messages
		for e = UGM_getEdges(n,edgeStruct)
			if n == edgeEnds(e,1)
				tmp = (edgePot(:,:,e) .* omsg(:,:,e+nEdges)') .^ (1 / c_ia(e));
			else
				tmp = (edgePot(:,:,e)' .* omsg(:,:,e)') .^ (1 / c_ia(e+nEdges));
			end
			if maximize
				% Max-product
				m_ai = max(tmp,[],2);
			else
				% Sum-product
				m_ai = sum(tmp,2);
			end
			if n == edgeEnds(e,1)
				imsg(:,e) = m_ai;
			else
				imsg(:,e+nEdges) = m_ai;
			end
		end
		
		% Beliefs
		newb = ones(nStates,1);
		for e = UGM_getEdges(n,edgeStruct)
			if n == edgeEnds(e,1)
				newb = newb .* imsg(:,e) .^ (c_ia(e) / c_i(n));
			else
				newb = newb .* imsg(:,e+nEdges) .^ (c_ia(e+nEdges) / c_i(n));
			end
		end
		z = sum(newb);
		if z > 0
			newb = newb / z;
		end
		bel(:,n) = newb;
		
		% Outgoing messages
		for e = UGM_getEdges(n,edgeStruct)
			if n == edgeEnds(e,1)
				tmp1 = (edgePot(:,:,e) .* omsg(:,:,e+nEdges)') .^ (-auxCount(e) / c_ia(e));
				tmp2 = (bel(:,n) ./ imsg(:,e)) .^ c_a(e);
				omsg(:,:,e) = tmp1 .* repmat(tmp2,1,nStates);
			else
				tmp1 = (edgePot(:,:,e)' .* omsg(:,:,e)') .^ (-auxCount(e+nEdges) / c_ia(e+nEdges));
				tmp2 = (bel(:,n) ./ imsg(:,e+nEdges)) .^ c_a(e);
				omsg(:,:,e) = tmp1 .* repmat(tmp2,1,nStates);
			end
		end
	end
	
	% Check convergence
	if abs(bel(:)-bel_old(:)) < convTol
		break
	end
	
	% Update old beliefs
	bel_old = bel;

% 	% Temp messages
% 	tmp_i = zeros(maxState,nNodes);
% 	tmp_o = zeros(maxState,nEdges*2);
% 	for e = 1:nEdges
% 		n1 = edgeEnds(e,1);
% 		n2 = edgeEnds(e,2);
% 		
% 		% Incoming
% 		pot = edgePot(1:nStates(n1),1:nStates(n2),e);
% 		if maximize
% 			% Max-product
% 			tmp_i(1:nStates(n1),e) = max(pot.*repmat(msg_o(1:nStates(n2),e+nEdges)',nStates(n1),1), [], 2);
% 			tmp_i(1:nStates(n2),e+nEdges) = max(pot'.*repmat(msg_o(1:nStates(n1),e)',nStates(n2),1), [], 2);
% 		else
% 			% Sum-product
% 			tmp_i(1:nStates(n1),e) = pot * msg_o(1:nStates(n2),e+nEdges);
% 			tmp_i(1:nStates(n2),e+nEdges) = pot' * msg_o(1:nStates(n1),e);
% 		end
% 		
% 		% Outgoing
% 		%  n1
% 		tmp = nodePot(n1,1:nStates(n1))';
% 		for e1 = UGM_getEdges(n1,edgeStruct)
% 			if e1 ~= e
% 				if n1 == edgeEnds(e1,1)
% 					tmp = tmp .* msg_i(1:nStates(n1),e1);
% 				else
% 					tmp = tmp .* msg_i(1:nStates(n1),e1+nEdges);
% 				end
% 			end
% 		end
% 		tmp_o(1:nStates(n1),e) = tmp;
% 		%  n2
% 		tmp = nodePot(n2,1:nStates(n2))';
% 		for e2 = UGM_getEdges(n2,edgeStruct)
% 			if e2 ~= e
% 				if n2 == edgeEnds(e2,1)
% 					tmp = tmp .* msg_i(1:nStates(n2),e2);
% 				else
% 					tmp = tmp .* msg_i(1:nStates(n2),e2+nEdges);
% 				end
% 			end
% 		end
% 		tmp_o(1:nStates(n2),e+nEdges) = tmp;
% 	end
% 	
% 	% New messages
% 	for e = 1:nEdges
% 		n1 = edgeEnds(e,1);
% 		n2 = edgeEnds(e,2);
% 		
% 		% exponent variables
% 		q1 = (1 - nodeCount(n1)) / length(UGM_getEdges(n1,edgeStruct));
% 		q2 = (1 - nodeCount(n2)) / length(UGM_getEdges(n2,edgeStruct));
% 		d1 = edgeCount(e) - q1 + 1;
% 		d2 = edgeCount(e) - q2 + 1;
% 		
% 		% Incoming
% 		%  n1
% 		newm = tmp_i(1:nStates(n1),e).^(edgeCount(e)/d1) .* ...
% 			   tmp_o(1:nStates(n1),e).^((q1-edgeCount(e))/d1);
% 		newm(~isfinite(newm)) = 0;
% 		z = sum(newm);
% 		if z > 0
% 			newm = newm ./ z;
% 		end
% 		msg_i(1:nStates(n1),e) = newm;
% 		%  n2
% 		newm = tmp_i(1:nStates(n2),e+nEdges).^(edgeCount(e)/d2) .* ...
% 			   tmp_o(1:nStates(n2),e+nEdges).^((q2-edgeCount(e))/d2);
% 		newm(~isfinite(newm)) = 0;
% 		z = sum(newm);
% 		if z > 0
% 			newm = newm ./ z;
% 		end
% 		msg_i(1:nStates(n2),e+nEdges) = newm;
% 
% 		% Outgoing
% 		%  n1
% 		newm = tmp_i(1:nStates(n1),e).^((q1-1)/d1) .* ...
% 			   tmp_o(1:nStates(n1),e).^(1/d1);
% 		newm(~isfinite(newm)) = 0;
% 		z = sum(newm);
% 		if z > 0
% 			newm = newm ./ z;
% 		end
% 		msg_o(1:nStates(n1),e) = newm;
% 		%  n2
% 		newm = tmp_i(1:nStates(n2),e+nEdges).^((q2-1)/d2) .* ...
% 			   tmp_o(1:nStates(n2),e+nEdges).^(1/d2);
% 		newm(~isfinite(newm)) = 0;
% 		z = sum(newm);
% 		if z > 0
% 			newm = newm ./ z;
% 		end
% 		msg_o(1:nStates(n2),e+nEdges) = newm;
% 	end
	
% 	% Check for NaNs
% 	[min(msg_i(:)) max(msg_i(:)) min(msg_o(:)) max(msg_o(:))]
% 	if any(isnan(msg_i(:))) || any(isnan(msg_o(:)))
% 		error('Found NaN values\n')
% 	end

% 	% Damping
% 	if momentum < 1
% 		msg_i = old_msg_i.^(1-momentum) .* msg_i.^momentum;
% 		msg_o = old_msg_o.^(1-momentum) .* msg_o.^momentum;
% 	end
	
% 	% Check convergence
% 	if (abs(imsg(:)-imsg_old(:)) + abs(omsg(:)-omsg_old(:))) < convTol
% 		break
% 	end
% 	
% 	% Swap old/new outgoing messages
% 	imsg_old = imsg;
% 	omsg_old = omsg;
end
if t == edgeStruct.maxIter
	fprintf('CountBP did not converge after %d iterations\n',edgeStruct.maxIter);
end
fprintf('CountBP stopped after %d iterations\n',t);
