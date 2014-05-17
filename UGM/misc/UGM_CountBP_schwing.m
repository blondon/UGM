function [imsg,omsg] = UGM_CountBP_schwing(nodePot,edgePot,nodeCount,edgeCount,edgeStruct,momentum,convTol,maximize)

[nNodes,nState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;

% For simplicity, all variables must have same number of states.
assert(all(edgeStruct.nStates == nState), 'UGM_CountBP: All variables must have same number of states.')

nodePot = log(nodePot);
edgePot = log(edgePot);

% Init messages
imsg = zeros(nState,nEdges*2); % incoming: e -> n
omsg = zeros(nState,nEdges*2); % outgoing: n -> e
itmp = zeros(nState,nEdges*2); % incoming: e -> n
otmp = zeros(nState,nEdges*2); % outgoing: n -> e
imsg_old = imsg;
omsg_old = omsg;

% Main loop
for i = 1:edgeStruct.maxIter
    
    % Iterate over nodes, in sequence
    for n = 1:nNodes
        
        % Neighbors of n
        neighbs = UGM_getEdges(n,edgeStruct);
        
        hatNodeCount = nodeCount(n) + sum(edgeCount(neighbs));

        
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

            iexp = 1/edgeCount(e);
            if maximize
                itmp(:,eidx) = max(iexp*(pot_ij + repmat(msg',nState,1)), [], 2);
            else
                itmp(:,eidx) = logsumexp(iexp*bsxfun(@plus, pot_ij, msg')')';
            end
            % Outgoing
            prod = nodePot(n,:)';
            for e_ = neighbs
                if n == edgeEnds(e_,1)
                    prod = prod + imsg(:,e_);
                else
                    prod = prod + imsg(:,e_+nEdges);
                end
            end
            otmp(:,eidx) = prod * (edgeCount(e)/hatNodeCount) - imsg(:, eidx);
        end
        
        % Update messages
        for e = neighbs
            if n == edgeEnds(e,1)
                eidx = e;
            else
                eidx = e + nEdges;
            end
            imsg(:,eidx) = logNormalizeMessage(itmp(:,eidx) * edgeCount(e));
            omsg(:,eidx) = logNormalizeMessage(otmp(:,eidx));
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
