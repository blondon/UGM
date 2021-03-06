function [new_msg] = UGM_TRBP(nodePot,edgePot,edgeStruct,maximize,mu)

[nNodes,maxState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = double(edgeStruct.nStates);

% Initialize
old_msg = zeros(maxState,nEdges*2);
new_msg = zeros(maxState,nEdges*2);
for e = 1:nEdges
	n1 = edgeEnds(e,1);
	n2 = edgeEnds(e,2);
	new_msg(1:nStates(n2),e) = 1/nStates(n2); % Message from n1 => n2
	new_msg(1:nStates(n1),e+nEdges) = 1/nStates(n1); % Message from n2 => n1
end

for i = 1:edgeStruct.maxIter
	for n = 1:nNodes
		% Find Neighbors
		edges = E(V(n):V(n+1)-1);
		
		% Send a message to each neighbor
		for e = edges(:)'
			n1 = edgeEnds(e,1);
			n2 = edgeEnds(e,2);
			
			if n == edgeEnds(e,2)
				pot_ij = edgePot(1:nStates(n1),1:nStates(n2),e);
			else
				pot_ij = edgePot(1:nStates(n1),1:nStates(n2),e)';
			end
			
			% Adjust edge potential by edge appearnce probability
			pot_ij = pot_ij.^(1/mu(e));
			
			% Compute temp = product of all incoming msgs except j
			%   to the power of the edge appearance probability,
			%   divided by msg from j to the (1 - edge appearnce prob)
			temp = nodePot(n,1:nStates(n))';
			for e2 = edges(:)'
				if n == edgeEnds(e2,2)
					incoming = new_msg(1:nStates(n),e2);
				else
					incoming = new_msg(1:nStates(n),e2+nEdges);
				end
				if e ~= e2
					temp = temp .* incoming.^mu(e2);
				else
					temp = temp ./ incoming.^(1-mu(e2));
				end
			end
			
			% Compute new message
			if maximize
				% These functions don't exist!
				if edgeStruct.useMex
					newm = max_mult(pot_ij,temp);
				else
					newm = max_multM(pot_ij,temp);
				end
			else
				newm = pot_ij * temp;
			end
			% Might get NaN or Inf values
			newm(~isfinite(newm)) = 0;
			
			% Safe normalize
			if sum(newm) > 0
				newm = newm ./ sum(newm);
			end
			
			if n == edgeEnds(e,2);
				new_msg(1:nStates(n1),e+nEdges) = newm;
			else
				new_msg(1:nStates(n2),e) = newm;
			end
		end
	end
	
% 	% Check for instability
% 	if any(isnan(new_msg(:))) || any(isinf(new_msg(:)))
% 		break
% 	end
	
	%fprintf('%d, %.4e\n',i,sum(abs(new_msg(:)-old_msg(:))));
	% Check convergence
	if sum(abs(new_msg(:)-old_msg(:))) < 1e-4
		break;
	end
	
	old_msg = new_msg;
end
if i == edgeStruct.maxIter
	fprintf('TRBP did not converge\n');
end
% fprintf('TRBP stopped after %d iterations\n',i);

