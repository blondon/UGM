function [nodeCount,edgeCount,auxCount] = UGM_ConvexBetheCounts(edgeStruct,kappa,tgt,verbose)
%
% Computes the counting numbers for the Bethe approximation.
%
% edgeStruct : edge structure
% kappa : desired modulus of convexity (def: 0)
% tgt : target counting numbers:
%		1 = Bethe (def) (Meshi et al., UAI'09)
%		2 = Uniform c_e=1 (Hazan & Shashua, UAI'08)
%		3 = TRW
% verbose : display log? (def: 1)
%
% nodeCount : (nNodes x 1) vector of node counting numbers
% edgeCount : (nEdges x 1) vector of edge counting numbers
% auxCount : (2*nEdges x 1) vector of auxiliary counting numbers

if ~exist('kappa','var') || isempty(kappa)
	kappa = 0;
end
if ~exist('tgt','var') || isempty(tgt)
	tgt = 1;
end
if ~exist('verbose','var') || isempty(verbose)
	verbose = 1;
end

tolCon = 1e-8;
tolValid = 1e-8;

% Compute target counting numbers
if tgt == 1
	[tgtNode,tgtEdge] = UGM_BetheCounts(edgeStruct);
elseif tgt == 2
	tgtNode = [];
	tgtEdge = ones(edgeStruct.nEdges,1);
elseif tgt == 3
	[tgtNode,tgtEdge] = UGM_TRBPCounts(edgeStruct);
else
	error('tgt must be 1 (Bethe), 2 (unif c_e) or 3 (TRW)')
end

% Try non-slack QP
[nodeCount,edgeCount,auxCount,~,exitflag] = solveQP(edgeStruct,kappa,tgtNode,tgtEdge,tolCon);
slack = [];

% Check feasibility
if exitflag == -2
	if verbose
		fprintf('QP is infeasible. Trying slackened version ...\n');
	end
	% Try slackened QP
	[nodeCount,edgeCount,auxCount,slack,~,exitflag] = solveSlackQP(edgeStruct,kappa,tgtNode,tgtEdge,tolCon);
	if exitflag == -2
		error('Slackened QP is infeasible. Try a different value for kappa.\n');
		return;
	end
end

if verbose
	
	% Convexity
	nEdges = edgeStruct.nEdges;
	min_a_e = inf;
	for e = 1:nEdges
		v = edgeCount(e);
		v = v - auxCount(e);
		v = v - auxCount(e+nEdges);
		if min_a_e > v
			min_a_e = v;
		end
	end
	fprintf('Solution is at least (%f-strongly) convex\n', min_a_e/3);

	% L2 distance^2 to target counts
	if tgt == 2
		fprintf('MSE target counts: %f\n', mean((tgtEdge-edgeCount).^2));
		fprintf('Max target counts: %f\n', max(abs(tgtEdge-edgeCount)));
	else
		fprintf('MSE target counts: %f\n', mean(([tgtNode;tgtEdge]-[nodeCount;edgeCount]).^2));
		fprintf('Max target counts: %f\n', max(abs([tgtNode;tgtEdge]-[nodeCount;edgeCount])));
	end
	
	% L2 distance^2 to variable-validity
	if ~isempty(slack)
		fprintf('MSE variable-validity: %f\n', mean(slack.^2));
		fprintf('Max variable-validity: %f\n', max(abs(slack)));
	end

	% LB constraints
	if any(edgeCount < 0) || any(auxCount < 0)
		fprintf('LB constraints not satisfied\n');
	end

end

%% Non-slackened QP
function [nodeCount,edgeCount,auxCount,fval,exitflag] = solveQP(edgeStruct,kappa,tgtNode,tgtEdge,tolCon)

% Dimensions
nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
nCnt = nNodes + nEdges;
nAux = 2*nEdges;
nVar = nCnt + nAux;

% Setup QP objective
if isempty(tgtNode)
	H = sparse(nNodes+1:nCnt,nNodes+1:nCnt,ones(nEdges,1),nVar,nVar);
	f = [zeros(nNodes,1) ; -2*tgtEdge ; zeros(nAux,1)];
elseif isempty(tgtEdge)
	H = sparse(1:nNodes,1:nNodes,ones(nNodes,1),nVar,nVar);
	f = [-2*tgtNode ; zeros(nEdges+nAux,1)];
else
	H = sparse(1:nCnt,1:nCnt,ones(nCnt,1),nVar,nVar);
	f = [-2*tgtNode ; -2*tgtEdge ; zeros(nAux,1)];
end

% Setup constraints
I = zeros(3*nEdges,1);
J = zeros(3*nEdges,1);
V = zeros(3*nEdges,1);
c = 0; i = 0;
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
A = sparse(I,J,V,c,nVar);
b = -ones(c,1) * kappa * 3;

[Aeq,beq] = variableValid(edgeStruct,nVar);

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


%% QP with count + slack variable-validity
function [nodeCount,edgeCount,auxCount,slack,fval,exitflag] = solveSlackQP(edgeStruct,kappa,tgtNode,tgtEdge,tolCon)

% Dimensions
nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
nCnt = nNodes + nEdges;
nAux = 2*nEdges;
nSlack = nNodes;
nVar = nCnt + nAux + nSlack;

% Setup QP objective
if isempty(tgtNode)
	H = sparse([nNodes+1:nCnt,nCnt+nAux+1:nVar],[nNodes+1:nCnt,nCnt+nAux+1:nVar],...
				ones(nEdges+nSlack,1),nVar,nVar);
	f = [zeros(nNodes,1) ; -2*tgtEdge ; zeros(nAux+nSlack,1)];
elseif isempty(tgtEdge)
	H = sparse([1:nNodes,nCnt+nAux+1:nVar],[1:nNodes,nCnt+nAux+1:nVar],...
				ones(nNodes+nSlack,1),nVar,nVar);
	H = sparse(1:nNodes,1:nNodes,ones(nNodes,1),nVar,nVar);
	f = [-2*tgtNode ; zeros(nEdges+nAux+nSlack,1)];
else
	H = sparse([1:nCnt,nCnt+nAux+1:nVar],[1:nCnt,nCnt+nAux+1:nVar],...
				ones(nCnt+nSlack,1),nVar,nVar);
	f = [-2*tgtNode ; -2*tgtEdge ; zeros(nAux+nSlack,1)];
end

% Setup constraints
I = zeros(3*nEdges,1);
J = zeros(3*nEdges,1);
V = zeros(3*nEdges,1);
c = 0; i = 0;
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
A = sparse(I,J,V,c,nVar);
b = -ones(c,1) * kappa * 3;

% Variable-validity slack constraints
I = zeros(2*nNodes+2*nEdges,1);
J = zeros(2*nNodes+2*nEdges,1);
V = zeros(2*nNodes+2*nEdges,1);
i = 0;
% foall v, z_v + c_v + sum_{e : v in e} c_e = 1
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
	i = i + 1;
	I(i) = n;
	J(i) = nCnt + nAux + n;
	V(i) = 1;
end
Aeq = sparse(I,J,V,nNodes,nVar);
beq = ones(nNodes,1);

lb = [-inf(nCnt,1) ; zeros(nAux,1) ; -inf(nSlack,1)];
ub = [inf(nCnt,1) ; inf(nAux,1) ; inf(nSlack,1)];

% Solve QP
options = optimset('Algorithm','interior-point-convex',...
				   'Display','off','TolCon',tolCon);
[x,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);

% Output
nodeCount = x(1:nNodes);
edgeCount = x(nNodes+1:nCnt);
auxCount = x(nCnt+1:nCnt+nAux);
slack = x(nCnt+nAux+1:end);


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

