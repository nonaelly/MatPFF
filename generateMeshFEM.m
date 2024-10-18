function [node, elem, nodeBou, elemBou] = generateMeshFEM(geoType, varargin)
switch geoType
    case 'cyl'
        R1 = varargin{1};
        R2 = varargin{2};
        O = varargin{3};
        numR = varargin{4};
        numTheta = varargin{5};
        dim = 2;
        node = zeros((numR+1) * (numTheta+1), 3);
        elem = zeros(numR * numTheta, 6);
        elemBou = cell(4, 1);

        ind = 1;
        for i = 1:numTheta+1
            for j = 1:numR+1
                R = R2 - (R2-R1) * (j - 1) / numR;
                theta = 90 - (90-0) * (i - 1) / numTheta;
                x = R * cosd(theta);
                y = R * sind(theta);
                %                 node((i-1)*(numR+1) + j, :) = [x, y];
                node(ind, :) = [ind, x, y];
                ind = ind + 1;
            end
        end
        node(2:3) = node(2:3) + O;
        nodeBou = cell(4, 1);
        nodeBou{1} = node(numR+1 : numR+1 : end, :);
        nodeBou{2} = node((numR+1)*numTheta+1 : end, :);
        nodeBou{3} = node(1 : numR+1 : end, :);
        nodeBou{4} = node(1:numR+1, :);

        ind = 1;
        for i = 1:numTheta
            for j = 1:numR
                ind1 = j + (i-1)*(numR+1);
                ind2 = ind1 + 1;
                ind3 = j + i*(numR+1) + 1;
                ind4 = ind3 - 1;
                elem(ind, :) = [ind, 1, ind1, ind2, ind3, ind4];
                if j == 1 % line 3
                    elemBou{3} = [elemBou{3}; ind, 4];
                end
                if j == numR % line 1
                    elemBou{1} = [elemBou{1}; ind, 2];
                end
                if i == 1 % line 4
                    elemBou{4} = [elemBou{4}; ind, 1];
                end
                if i == numTheta % line 2
                    elemBou{2} = [elemBou{2}; ind, 3];
                end
                ind = ind + 1;
            end
        end
        figure;
        hold on;
        for i = 1:size(elem, 1)
            nodes = node(elem(i, 3:end), 2:3);
            fill(nodes(:, 1), nodes(:, 2), 'w', 'EdgeColor', 'k');
        end
        plot(node(:, 2), node(:, 3), 'ro', 'MarkerFaceColor', 'r');
        title('Mesh Plot');
        xlabel('X');
        ylabel('Y');
        axis equal;
    case 'rect'
        width = varargin{1};
        height = varargin{2};
        O = varargin{3};
        numX = varargin{4};
        numY = varargin{5};

        % 生成节点
        node = zeros((numX+1) * (numY+1), 3);
        ind = 1;
        for i = 1:numY+1
            for j = 1:numX+1
                x = width * (j-1) / numX;
                y = height * (i-1) / numY;
                node(ind, :) = [ind, x, y];
                ind = ind + 1;
            end
        end
        node(2:3) = node(2:3) + O;

        elem = zeros(numX * numY, 6);
        elemBou = cell(4, 1);
        ind = 1;
        for i = 1:numY
            for j = 1:numX
                ind1 = j + (i-1)*(numX+1);
                ind2 = ind1 + 1;
                ind3 = ind1 + (numX+1);
                ind4 = ind3 + 1;

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
        nodeBou{1} = node(1:numX+1, :);
        nodeBou{2} = node((numX+1):numX+1:end, :);
        nodeBou{3} = node(end-numX:end, :);
        nodeBou{4} = node(1:numX+1:end, :);

        % 绘制网格
        figure;
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