function [node, elem, nodeBou, elemBou] = generateMeshFEM_CZM(geoType, varargin)
switch geoType
    case 'cyl'

    case 'rect'
        width = varargin{1};
        height = varargin{2};
        O = varargin{3};
        numX = varargin{4};
        numY = varargin{5};

        % nodes in CZM are double
        node = zeros((numX+2) * (numY+1), 3);
        ind = 1;
        for i = 1:numY+1
            for j = 1:numX+2
                if j <= numX/2+1
                    x = width * (j-1) / numX;
                else
                    x = width * (j-2) / numX;
                end
                y = height * (i-1) / numY;
                node(ind, :) = [ind, x, y];
                ind = ind + 1;
            end
        end
        node(:, 2:3) = node(:, 2:3) + O;

        elem = zeros(numX * numY, 6);
        elemBou = cell(4, 1);
        ind = 1;
        for i = 1:numY
            for j = 1:numX
                if j <= numX/2
                    ind1 = j + (i-1)*(numX+2);
                    ind2 = ind1 + 1;
                    ind3 = ind1 + (numX+2);
                    ind4 = ind3 + 1;
                else
                    ind1 = j + 1 + (i-1)*(numX+2);
                    ind2 = ind1 + 1;
                    ind3 = ind1 + (numX+2);
                    ind4 = ind3 + 1;
                end


                if i == 1
                    elemBou{1} = [elemBou{1}; ind, 1];
                end
                if j == numX
                    elemBou{2} = [elemBou{2}; ind, 2];
                end
                if i == numX
                    elemBou{3} = [elemBou{3}; ind, 3];
                end
                if j == 1 %
                    elemBou{4} = [elemBou{4}; ind, 4];
                end

                elem(ind, :) = [ind, 1, ind1, ind2, ind4, ind3];
                ind = ind + 1;
            end
        end

        nodeBou = cell(4, 1);
        nodeBou{1} = node(1:numX+2, :);
        nodeBou{2} = node((numX+2):numX+2:end, :);
        nodeBou{3} = node(end-numX-1:end, :);
        nodeBou{4} = node(1:numX+2:end, :);

        % 绘制网格
        figure(1);
        hold on;
        for i = 1:size(elem, 1)
            nodes = node(elem(i, 3:end), 2:3);
            fill(nodes(:, 1), nodes(:, 2), 'w', 'EdgeColor', 'k');
        end
        plot(node(:, 2), node(:, 3), 'ro', 'MarkerFaceColor', 'r');
        title('Rectangular Mesh Plot');
        xlabel('X');
        ylabel('Y');
        axis equal;
end
end