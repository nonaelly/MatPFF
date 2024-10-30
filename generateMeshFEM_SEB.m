function [node, elem, nodeBou, elemBou] = generateMeshFEM_SEB(dx, numY)
w = 100;
a = 19;
t = 75;
b1 = 162;
b2 = 26;

% Cohesive node set
CNode = [];

% ---------------------------------------
% Generate nodes
node = [];
nodeBou = cell(4, 1);

% Entity node
x0 = 0;
x = x0;
ind = 1;
indX = 0;
% (0,0) - (162,0)
while x0 <= b1
    for j = 1:numY(1)
        y = a * (j-1) / numY(1);
        [node] = [node; ind, x, y];
        ind = ind + 1;
    end
    for j = 1 : numY(2)+1
        y = (w-a) * (j-1) / numY(2) + a;
        [node] = [node; ind, x, y];
        if x0 == 0
            [CNode] = [CNode; ind, x, y];
        end
        ind = ind + 1;
    end

    x = x0 + dx;
    x0 = x;

    if x0 > b1 && x0 - dx < b1
        x = b1;
        x0 = x;
    end

    indX = indX + 1;
    if mod(indX, 5) == 0 && dx < 3
        dx = dx * 2;
    end

end
% (162,0) - (188,0)
for i = 1:8
    x = b1 + b2 * i / 8;
    for j = 1:numY(1)+1
        y = a * (j-1) / numY(1);
        [node] = [node; ind, x, y];
        ind = ind + 1;
    end
    for j = 1 : numY(2)
        y = (w-a) * j / numY(2) + a;
        [node] = [node; ind, x, y];
        ind = ind + 1;
    end
    indX = indX + 1;
end

nodeBou{1} = node(1:sum(numY)+1:end, :);
nodeBou{2} = node(end-sum(numY):end, :);
nodeBou{3} = node(sum(numY)+1:sum(numY)+1:end, :);

% Cohesive node
% (0,0) - (0,19)
for j = 1 : numY(2) + 1
    x = 0;
    y = (w-a) * (j-1) / numY(2) + a;
    [node] = [node; ind, x, y];
    [CNode] = [CNode; ind, x, y];
    nodeBou{4} = [nodeBou{4}; ind, x, y];
    ind = ind + 1;
end

% ---------------------------------------
% Generate elements
numX = indX - 1;
elem = [];
elemBou = cell(4, 1);

% Entity elements
ind = 1;
for i = 1:numX
    for j = 1:sum(numY)
        ind1 = j + (i-1)*(sum(numY)+1);
        ind2 = ind1 + (sum(numY)+1);
        ind3 = ind2 + 1;
        ind4 = ind1 + 1;

        if j == 1
            elemBou{1} = [elemBou{1}; ind, 1];
        end
        if i == numX
            elemBou{2} = [elemBou{2}; ind, 2];
        end
        if j == sum(numY)
            elemBou{3} = [elemBou{3}; ind, 3];
        end

        % The fourth boundary is not here.

        % The second index: element type.
        [elem] = [elem; ind, 1, ind2, ind3, ind4, ind1];
        ind = ind + 1;
    end
end

% Cohesive elements
for j = 1:numY(2)
    ind1 = j;
    ind2 = ind1 + (numY(2) + 1);
    ind3 = ind2 + 1;
    ind4 = ind1 + 1;

    elemBou{4} = [elemBou{4}; ind, 4];

    % The second index: element type.
    [elem] = [elem; ind, 2, CNode(ind2, 1), CNode(ind3, 1), CNode(ind4, 1), CNode(ind1, 1)];
    ind = ind + 1;
end

% figure(1);
% hold on;
% for i = 1:size(elem, 1)
%     nodes = node(elem(i, 3:end), 2:3);
%     fill(nodes(:, 1), nodes(:, 2), 'w', 'EdgeColor', 'k');
%     if elem(i, 2) == 2 % CZM
%         plot(node(elem(i, 3:end), 2), node(elem(i, 3:end), 3), 'ro', 'MarkerFaceColor', 'r');
%     end
% end
% 
% title('single-edge notched beam (SE(B)) test');
% xlabel('X');
% ylabel('Y');
% axis equal;

end