function [fixNode, nodeForce] = generateBC_CZM_Bilinear(nodeBou, node, ubar)

fixNode = [];
nodeForce = [];

% (162, 0)
m = find(ismember(node{1}, [162, 0], 'rows'));
% [fixNode] = [fixNode; m, 1, 0; m, 2, 0];
[fixNode] = [fixNode; m, 2, 0];

% (0, 100) P
m = find(ismember(node{1}, [0, 100], 'rows'));
[fixNode] = [fixNode; m, 2, ubar; size(node{1}, 1) + size(node{2}, 1)/2, 2, ubar];

% part2 line4 X
line = nodeBou{2}{4} + size(node{1}, 1);
[fixNode] = [fixNode; line(:, 1), 1*ones(size(line, 1), 1), zeros(size(line, 1), 1)];

fixNode = sortrows(fixNode, 1);
end