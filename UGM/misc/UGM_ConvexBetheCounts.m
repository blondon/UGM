function [nodeCount,edgeCount,auxCount] = UGM_ConvexBetheCounts(edgeStruct,kappa,minKappa)
%
% Computes the counting numbers for the Bethe approximation.
%
% edgeStruct : edge structure
% kappa : minimum modulus of convexity (def: 0)
%
% nodeCount : (nNodes x 1) vector of node counting numbers
% edgeCount : (nEdges x 1) vector of edge counting numbers
% auxCount : (2*nEdges x 1) vector of auxiliary counting numbers

if nargin < 2
	kappa = 0;
end
if nargin < 3
	minKappa = 1e-2;
end

tolCon = 1e-8;
tolValid = 1e-8;
tolConvex = 1e-8;

% Compute Bethe counting numbers
[nodeCount_b,edgeCount_b] = UGM_BetheCounts(edgeStruct);
betheCount = [nodeCount_b ; edgeCount_b];

% % TEMP: for TRBP testing
% nodeCount_b = ones(edgeStruct.nNodes,1);
% edgeCount_b = edgeStruct.edgeDist;
% for e = 1:edgeStruct.nEdges
% 	n1 = edgeStruct.edgeEnds(e,1);
% 	n2 = edgeStruct.edgeEnds(e,2);
% 	nodeCount_b(n1) = nodeCount_b(n1) - edgeStruct.edgeDist(e);
% 	nodeCount_b(n2) = nodeCount_b(n2) - edgeStruct.edgeDist(e);
% end
% betheCount = [nodeCount_b ; edgeCount_b];


% Try non-slack QP
[nodeCount,edgeCount,auxCount,~,exitflag] = solveQP(edgeStruct,kappa,betheCount,tolCon);
slack = [];

% Check feasibility
if exitflag == -2
	fprintf('QP is infeasible. Trying slackened version with minKappa=%f ...\n',minKappa);
	
	% Try slackened QP
	[nodeCount,edgeCount,auxCount,slack,~,exitflag] = solveSlackQP(edgeStruct,kappa,minKappa,betheCount,tolCon);
	if exitflag == -2
		fprintf('Slackened QP is infeasible. Try a different value for kappa,minKappa.\n');
	end
end

% Objective value
fprintf('Fit: %f\n', norm([nodeCount;edgeCount]-betheCount,2)^2);

% Verify
verifySolution(edgeStruct,kappa,nodeCount,edgeCount,auxCount,slack,tolValid);


%% Non-slackened QP
function [nodeCount,edgeCount,auxCount,fval,exitflag] = solveQP(edgeStruct,kappa,betheCount,tolCon)

% Dimensions
nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
nCnt = nNodes + nEdges;
nAux = 2*nEdges;
nVar = nCnt + nAux;

% Setup QP variables
H = sparse(1:nCnt,1:nCnt,ones(nCnt,1),nVar,nVar);
f = [betheCount ; zeros(nAux,1)];

I = zeros(nNodes+5*nEdges,1);
J = zeros(nNodes+5*nEdges,1);
V = zeros(nNodes+5*nEdges,1);
c = 0; i = 0;
for n = 1:nNodes
	c = c + 1;
	i = i + 1;
	I(i) = c;
	J(i) = n;
	V(i) = -1;
	for e = UGM_getEdges(n,edgeStruct)
		i = i + 1;
		I(i) = c;
		if n == edgeStruct.edgeEnds(e,1)
			J(i) = nCnt + e;
		else
			J(i) = nCnt + e + nEdges;
		end
		V(i) = -1;
	end
end
for e = 1:nEdges
	c = c + 1;
	i = i + 1;
	I(i) = c;
	J(i) = nNodes + e;
	V(i) = -1;
	i = i + 1;
	I(i) = c;
	J(i) = nCnt + e;
	V(i) = 1;
	i = i + 1;
	I(i) = c;
	J(i) = nCnt + e + nEdges;
	V(i) = 1;
end
A = sparse(I,J,V,nCnt,nVar);
b = -ones(nCnt,1) * kappa;

I = zeros(nNodes+2*nEdges,1);
J = zeros(nNodes+2*nEdges,1);
V = zeros(nNodes+2*nEdges,1);
i = 0;
for n = 1:nNodes
	i = i + 1;
	I(i) = n;
	J(i) = n;
	V(i) = 1;
	for e = UGM_getEdges(n,edgeStruct)
		i = i + 1;
		I(i) = n;
		J(i) = nNodes + e;
		V(i) = 1;
	end
end
Aeq = sparse(I,J,V,nNodes,nVar);
beq = ones(nNodes,1);

lb = [-inf(nCnt,1) ; zeros(nAux,1)];
ub = [inf(nCnt,1) ; inf(nAux,1)];

% Solve QP
options = optimset('Algorithm','interior-point-convex',...
				   'Display','off','TolCon',tolCon);
[x,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);

% Output
nodeCount = x(1:nNodes);
edgeCount = x(nNodes+1:nCnt);
auxCount = x(nCnt+1:end);


%% Slackened QP
function [nodeCount,edgeCount,auxCount,slack,fval,exitflag] = solveSlackQP(edgeStruct,kappa,minKappa,betheCount,tolCon)

% Dimensions
nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
nCnt = nNodes + nEdges;
nAux = 2*nEdges;
nSlack = nCnt;
nVar = nCnt + nAux + nSlack;

% Setup QP variables
H = sparse(1:nCnt,1:nCnt,ones(nCnt,1),nVar,nVar);
f = [betheCount ; zeros(nAux,1) ; ones(nSlack,1)];

I = zeros(nCnt + 2*nAux + nSlack,1);
J = zeros(nCnt + 2*nAux + nSlack,1);
V = zeros(nCnt + 2*nAux + nSlack,1);
c = 0; i = 0;
for n = 1:nNodes
	c = c + 1;
	i = i + 1;
	I(i) = c;
	J(i) = n;
	V(i) = -1;
	for e = UGM_getEdges(n,edgeStruct)
		i = i + 1;
		I(i) = c;
		if n == edgeStruct.edgeEnds(e,1)
			J(i) = nCnt + e;
		else
			J(i) = nCnt + e + nEdges;
		end
		V(i) = -1;
	end
	i = i + 1;
	I(i) = c;
	J(i) = nCnt + nAux + n;
	V(i) = -1;
end
for e = 1:nEdges
	c = c + 1;
	i = i + 1;
	I(i) = c;
	J(i) = nNodes + e;
	V(i) = -1;
	i = i + 1;
	I(i) = c;
	J(i) = nCnt + e;
	V(i) = 1;
	i = i + 1;
	I(i) = c;
	J(i) = nCnt + e + nEdges;
	V(i) = 1;
	i = i + 1;
	I(i) = c;
	J(i) = nCnt + nAux + e;
	V(i) = -1;
end
A = sparse(I,J,V,nCnt,nVar);
b = -ones(nCnt,1) * kappa;

I = zeros(nNodes+2*nEdges,1);
J = zeros(nNodes+2*nEdges,1);
V = zeros(nNodes+2*nEdges,1);
i = 0;
for n = 1:nNodes
	i = i + 1;
	I(i) = n;
	J(i) = n;
	V(i) = 1;
	for e = UGM_getEdges(n,edgeStruct)
		i = i + 1;
		I(i) = n;
		J(i) = nNodes + e;
		V(i) = 1;
	end
end
Aeq = sparse(I,J,V,nNodes,nVar);
beq = ones(nNodes,1);

lb = [-inf(nCnt,1) ; zeros(nAux,1) ; zeros(nSlack,1)];
ub = [inf(nCnt,1) ; inf(nAux,1) ; (kappa-minKappa)*ones(nSlack,1)];

% Solve QP
options = optimset('Algorithm','interior-point-convex',...
				   'Display','off','TolCon',tolCon);
[x,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);

% Output
nodeCount = x(1:nNodes);
edgeCount = x(nNodes+1:nCnt);
auxCount = x(nCnt+1:nCnt+nAux);
slack = x(nCnt+nAux+1:end);


%% Verify solution
function verifySolution(edgeStruct,kappa,nodeCount,edgeCount,auxCount,slack,tolValid)

% LB constraints
if any(edgeCount < 0) || any(auxCount < 0)
	fprintf('LB constraints not satisfied\n');
else
	fprintf('LB constraints satisfied\n');
end

% Validity constraints
satisfied = 1;
maxViolation = 0;
for n = 1:edgeStruct.nNodes
	val = nodeCount(n);
	for e = UGM_getEdges(n,edgeStruct)
		val = val + edgeCount(e);
	end
	if abs(val - 1) > tolValid
		violation = abs(val - 1) - tolValid;
		%fprintf('Node %d count exceeded validity tolerance %f by %f\n',n,tolValid,violation);
		if maxViolation < violation
			maxViolation = violation;
		end
		satisfied = 0;
	end
end
if ~satisfied
	fprintf('Solution is not variable-valid; max violation: %f\n',maxViolation);
else
	fprintf('Solution is variable-valid\n');
end

% Convexity
convexity = kappa;
if ~isempty(slack)
	convexity = convexity - max(slack);
end
fprintf('Solution is (%f-strongly) convex\n',convexity);

% % Convexity constraints
% satisfied = 1;
% maxViolation = 0;
% for n = 1:nNodes
% 	val = nodeCount(n);
% 	for e = UGM_getEdges(n,edgeStruct)
% 		if n == edgeStruct.edgeEnds(e,1)
% 			val = val + auxCount(e);
% 		else
% 			val = val + auxCount(e+nEdges);
% 		end
% 	end
% 	if val < kappa-tolConvex
% 		violation = kappa - tolConvex - val;
% 		%fprintf('Node %d count exceeded convexity tolerance %f by %f\n',n,tolConvex,violation);
% 		if maxViolation < violation
% 			maxViolation = violation;
% 		end
% 		satisfied = 0;
% 	end
% end
% for e = 1:nEdges
% 	val = edgeCount(e) - auxCount(e) - auxCount(e+nEdges);
% 	if val < kappa-tolConvex
% 		violation = kappa - tolConvex - val;
% 		%fprintf('Edge %d count exceeded convexity tolerance %f by %f\n',e,tolConvex,violation);
% 		if maxViolation < violation
% 			maxViolation = violation;
% 		end
% 		satisfied = 0;
% 	end
% end
% if ~satisfied
% 	fprintf('Solution is not (%f-strongly) convex; max violation: %f\n',kappa,maxViolation);
% else
% 	fprintf('Solution is (%f-strongly) convex\n',kappa);
% end
