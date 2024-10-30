function [fixNode, nodeForce] = generateBC_CZM_BilinearV2(nodeBou, node, ubar)

fixNode = [];
nodeForce = [];

% (162, 0)
m = find(ismember(node, [162, 0], 'rows'));
% [fixNode] = [fixNode; m, 1, 0; m, 2, 0];
[fixNode] = [fixNode; m, 2, 0];

% (0, 100) P
m = find(ismember(node, [0, 100], 'rows'));
[fixNode] = [fixNode; m(1), 2, ubar; m(2), 2, ubar];

% part2 line4 X
line = nodeBou{4};
[fixNode] = [fixNode; line(:, 1), 1*ones(size(line, 1), 1), zeros(size(line, 1), 1)];

fixNode = sortrows(fixNode, 1);
end