function error = UGM_checkLocalMarginalConsistency(nodeBel, edgeBel, edgeStruct)

error = zeros(2*edgeStruct.nEdges,2);
for e = 1:edgeStruct.nEdges
    eBel = edgeBel(:,:,e);
    
    error(e,1) = sum(abs(sum(eBel')-nodeBel(edgeStruct.edgeEnds(e,1),:)));
    error(e,2) = sum(abs(sum(eBel)-nodeBel(edgeStruct.edgeEnds(e,2),:)));
end

