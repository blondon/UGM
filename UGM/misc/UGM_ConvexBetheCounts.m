function [nodeCount,edgeCount,auxCount] = UGM_ConvexBetheCounts(edgeStruct,kappa,tgt,verbose)
%
% Computes the counting numbers for the Bethe approximation.
%
% edgeStruct : edge structure
% kappa : desired modulus of convexity (def: 0)
% tgt : target counting numbers: 1 = Bethe (def), 2 = TRW
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
	[nodeCount_b,edgeCount_b] = UGM_BetheCounts(edgeStruct);
	tgtCount = [nodeCount_b ; edgeCount_b];
elseif tgt == 2
	[nodeCount_b,edgeCount_b] = UGM_TRBPCounts(edgeStruct);
	tgtCount = [nodeCount_b ; edgeCount_b];
else
	error('tgt must be 1 (Bethe) or 2 (TRW)')
end

% Try non-slack QP
[nodeCount,edgeCount,auxCount,~,exitflag] = solveQP(edgeStruct,kappa,tgtCount,tolCon);
slack = [];

% Check feasibility
if exitflag == -2
	if verbose
		fprintf('QP is infeasible. Trying slackened version ...\n');
	end
	% Try slackened QP
	[nodeCount,edgeCount,auxCount,slack,~,exitflag] = solveSlackQP(edgeStruct,kappa,tgtCount,tolCon);
	if exitflag == -2
		error('Slackened QP is infeasible. Try a different value for kappa.\n');
		return;
	end
end

if verbose
	
	% Convexity
	nNodes = edgeStruct.nNodes;
	nEdges = edgeStruct.nEdges;
	minKappa = inf;
	for n = 1:nNodes
		v = nodeCount(n);
		for e = UGM_getEdges(n,edgeStruct)
			if n == edgeStruct.edgeEnds(e,1)
				v = v + auxCount(e);
			else
				v = v + auxCount(e+nEdges);
			end
		end
		if minKappa > v
			minKappa = v;
		end
	end
	for e = 1:nEdges
		v = edgeCount(e);
		v = v + auxCount(e);
		v = v + auxCount(e+nEdges);
		if minKappa > v
			minKappa = v;
		end
	end
	fprintf('Solution is (%f-strongly) convex\n',minKappa);

	% L2 distance^2 to target counts
	fprintf('MSE target counts: %f\n', mean((tgtCount-[nodeCount;edgeCount]).^2));
	fprintf('Max target counts: %f\n', max(abs(tgtCount-[nodeCount;edgeCount])));

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
function [nodeCount,edgeCount,auxCount,fval,exitflag] = solveQP(edgeStruct,kappa,tgtCount,tolCon)

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
function [nodeCount,edgeCount,auxCount,slack,fval,exitflag] = solveSlackQP(edgeStruct,kappa,tgtCount,tolCon)

% Dimensions
nNodes = double(edgeStruct.nNodes);
nEdges = double(edgeStruct.nEdges);
nCnt = nNodes + nEdges;
nAux = 2*nEdges;
nSlack = nNodes;
nVar = nCnt + nAux + nSlack;

% Setup QP variables
H = sparse([1:nCnt,nCnt+nAux+1:nVar],[1:nCnt,nCnt+nAux+1:nVar],...
			ones(nCnt+nSlack,1),nVar,nVar);
f = [-2*tgtCount ; zeros(nAux+nSlack,1)];

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

