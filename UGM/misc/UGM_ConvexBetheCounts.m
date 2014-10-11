function [nodeCount,edgeCount,auxCount] = UGM_ConvexBetheCounts(edgeStruct,kappa,minKappa,validity,tgt)
%
% Computes the counting numbers for the Bethe approximation.
%
% edgeStruct : edge structure
% kappa : desired modulus of convexity (def: 0)
% minKappa : minimum modulus of convexity, if first QP fails (def: 0.01)
% validity : 0 = no validity constraints (def), 1 = variable-valid, 2 = factor-valid
% tgt : target counting numbers: 1 = Bethe (def), 2 = TRW
%
% nodeCount : (nNodes x 1) vector of node counting numbers
% edgeCount : (nEdges x 1) vector of edge counting numbers
% auxCount : (2*nEdges x 1) vector of auxiliary counting numbers

if ~exist('kappa','var') || isempty(kappa)
	kappa = 0;
end
if ~exist('minKappa','var') || isempty(minKappa)
	minKappa = 1e-2;
end
if ~exist('validity','var') || isempty(validity)
	validity = 0;
end
if ~exist('tgt','var') || isempty(tgt)
	tgt = 1;
end

tolCon = 1e-8;
tolValid = 1e-8;

% Compute target counting numbers
if tgt == 1
	[nodeCount_b,edgeCount_b] = UGM_BetheCounts(edgeStruct);
	tgtCount = [nodeCount_b ; edgeCount_b];
elseif tgt == 2
	[nodeCount_b,edgeCount_b] = UGM_TRBPCounts(edgeStruct);
	tgtCount = [nodeCount_b ; edgeCount_b];
else
	error('tgt must be 1 (Bethe) or 2 (TRW)')
end

% Try non-slack QP
[nodeCount,edgeCount,auxCount,~,exitflag] = solveQP(edgeStruct,kappa,tgtCount,validity,tolCon);
slack = [];

% Check feasibility
if exitflag == -2
	fprintf('QP is infeasible. Trying slackened version with minKappa=%f ...\n',minKappa);
	
	% Try slackened QP
	[nodeCount,edgeCount,auxCount,slack,~,exitflag] = solveSlackQP(edgeStruct,kappa,minKappa,tgtCount,validity,tolCon);
	if exitflag == -2
		fprintf('Slackened QP is infeasible. Try a different value for kappa,minKappa.\n');
	end
end

% Objective value
fprintf('Fit: %f\n', norm([nodeCount;edgeCount]-tgtCount,2)^2);

% Verify
verifySolution(edgeStruct,kappa,nodeCount,edgeCount,auxCount,slack,tolValid);


%% Non-slackened QP
function [nodeCount,edgeCount,auxCount,fval,exitflag] = solveQP(edgeStruct,kappa,tgtCount,validity,tolCon)

% Dimensions
nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
nCnt = nNodes + nEdges;
nAux = 2*nEdges;
nVar = nCnt + nAux;

% Setup QP variables
H = sparse(1:nCnt,1:nCnt,ones(nCnt,1),nVar,nVar);
f = [-2*tgtCount ; zeros(nAux,1)];

I = zeros(nNodes+5*nEdges,1);
J = zeros(nNodes+5*nEdges,1);
V = zeros(nNodes+5*nEdges,1);
c = 0; i = 0;
% forall v, -(c_v + sum_{e : v in e} a_{v,e}) <= -k
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
% forall e, -(c_e - sum_{v : v in e} a_{v,e}) <= -k
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

switch(validity)
	case 1
		[Aeq,beq] = variableValid(edgeStruct,nVar);
	case 2
		[Aeq,beq] = factorValid(edgeStruct,nVar);
	otherwise
		Aeq = []; beq = [];
end

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
function [nodeCount,edgeCount,auxCount,slack,fval,exitflag] = solveSlackQP(edgeStruct,kappa,minKappa,tgtCount,validity,tolCon)

% Dimensions
nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
nCnt = nNodes + nEdges;
nAux = 2*nEdges;
nSlack = nCnt;
nVar = nCnt + nAux + nSlack;

% Setup QP variables
H = sparse(1:nCnt,1:nCnt,ones(nCnt,1),nVar,nVar);
f = [-2*tgtCount ; zeros(nAux,1) ; ones(nSlack,1)];

I = zeros(nCnt + 2*nAux + nSlack,1);
J = zeros(nCnt + 2*nAux + nSlack,1);
V = zeros(nCnt + 2*nAux + nSlack,1);
c = 0; i = 0;
% forall v, -(c_v + sum_{e : v in e} a_{v,e}) <= -(k - z_v)
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
% forall e, -(c_e - sum_{v : v in e} a_{v,e}) <= -(k - z_e)
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

switch(validity)
	case 1
		[Aeq,beq] = variableValid(edgeStruct,nVar);
	case 2
		[Aeq,beq] = factorValid(edgeStruct,nVar);
	otherwise
		Aeq = []; beq = [];
end

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

% Variable-validity
satisfied = 1;
maxViolation = 0;
for n = 1:edgeStruct.nNodes
	val = nodeCount(n);
	for e = UGM_getEdges(n,edgeStruct)
		val = val + edgeCount(e);
	end
	if abs(val-1) > tolValid
		violation = abs(val-1) - tolValid;
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

% Factor-validity
satisfied = 1;
maxViolation = 0;
for e = 1:edgeStruct.nEdges
	val = edgeCount(e);
	if abs(val-1) > tolValid
		violation = abs(val-1) - tolValid;
		if maxViolation < violation
			maxViolation = violation;
		end
		satisfied = 0;
	end
end
if ~satisfied
	fprintf('Solution is not factor-valid; max violation: %f\n',maxViolation);
else
	fprintf('Solution is factor-valid\n');
end

% Convexity
convexity = kappa;
if ~isempty(slack)
	convexity = convexity - max(slack);
end
fprintf('Solution is (%f-strongly) convex\n',convexity);


%% Variable-validity constraints
function [Aeq,beq] = variableValid(edgeStruct,nVar)

nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
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


%% Factor-validity constraints
function [Aeq,beq] = factorValid(edgeStruct,nVar)

nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
I = zeros(nEdges,1);
J = zeros(nEdges,1);
V = zeros(nEdges,1);
for e = 1:nEdges
	I(e) = e;
	J(e) = nNodes + e;
	V(e) = 1;
end
Aeq = sparse(I,J,V,nEdges,nVar);
beq = ones(nEdges,1);



