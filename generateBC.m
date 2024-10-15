function [fixNode, nodeForce] = generateBC(BC, nodeBou, elemBou, elem, node, pres)

fixNode = [];
nodeForce = [];
indP = 1;
for i = 1 : 4
    line = nodeBou{i};
    switch BC{i}
        case {'dy', 'y'}
            [fixNode] = [fixNode; line(:, 1), 2*ones(size(line, 1), 1), zeros(size(line, 1), 1)];
        case {'dx', 'x'}
            [fixNode] = [fixNode; line(:, 1), 1*ones(size(line, 1), 1), zeros(size(line, 1), 1)];
        case {'p', 'f'}
            if BC{i} == 'p'
                bElem = elemBou{i}(:, 1);
                face = elemBou{i}(:, 2);
                faceNode = [1,2;2,3;3,4;4,1];
                tempF = zeros(4, 3);
                for j = 1:size(bElem, 1)
                    bNode = elem(bElem(j), faceNode(face(j),:)); % 单元的边界节点
                    tempF(4*j-3,1) = bNode(1);
                    tempF(4*j-2,1) = bNode(1);
                    tempF(4*j-1,1) = bNode(2);
                    tempF(4*j-0,1) = bNode(2);
                    tempF(4*j-3,2) = 1;
                    tempF(4*j-2,2) = 2;
                    tempF(4*j-1,2) = 1;
                    tempF(4*j-0,2) = 2;
                    tempF(4*j-3,3) = 0.5*pres(indP)*(node(bNode(1),2)-node(bNode(2),2));
                    tempF(4*j-2,3) = 0.5*pres(indP)*(node(bNode(2),1)-node(bNode(1),1));
                    tempF(4*j-1,3) = 0.5*pres(indP)*(node(bNode(1),2)-node(bNode(2),2));
                    tempF(4*j-0,3) = 0.5*pres(indP)*(node(bNode(2),1)-node(bNode(1),1));
                end
                [nodeForce] = [nodeForce; tempF];
                indP = indP + 1;
            end
        case {'d'}
            for j = 1 : size(line, 1)
                [fixNode] = [fixNode; line(j, 1), 1, 0; line(j, 1), 2, 0];
            end
    end
end
fixNode = sortrows(fixNode, 1);
end