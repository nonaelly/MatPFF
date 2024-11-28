function FN = findFailureNodes(sVar, elemCoh, coNode)
FN = [];

% possible failure cohesive node: coNode
for i = 1 : size(coNode, 1)
    elemIdx = [];
    for j = 1 : size(elemCoh, 1)
        if ismember(coNode(i, 1), elemCoh(j, :))
            [elemIdx] = [elemIdx; j];
        end
    end

    failure = 1;
    for j = 1 : size(elemIdx, 1)
        if ismember(sVar.isFail{elemIdx(j)}, 0)
            failure = 0;
            break
        end
    end

    if failure
        [FN] = [FN; coNode(i, 1)];
    end
end
FN = sortrows(FN);
end