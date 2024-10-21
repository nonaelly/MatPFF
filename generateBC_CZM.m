function [fixNode, nodeForce] = generateBC_CZM(nodeBou, node, ubar)

fixNode = [];
nodeForce = [];

% (0,0)
[fixNode] = [fixNode; 1, 1, 0];

% (0,0.5)
% [fixNode] = [fixNode; 1 + size(node{1}, 1), 1, 0];

% part1 line1 Y
line = nodeBou{1}{1};
[fixNode] = [fixNode; line(:, 1), 2*ones(size(line, 1), 1), zeros(size(line, 1), 1)];

% part2 line3 Y
line = nodeBou{2}{3} + size(node{1}, 1);
[fixNode] = [fixNode; line(:, 1), 2*ones(size(line, 1), 1), ubar*ones(size(line, 1), 1)];

fixNode = sortrows(fixNode, 1);
end