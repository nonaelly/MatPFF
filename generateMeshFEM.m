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
        node(:, 2:3) = node(:, 2:3) + O;

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
    case 'SEB'
        w = 100;
        a = 19;
        t = 75;
        b1 = 162;
        b2 = 26;
        dx = varargin{1};
        numY = varargin{2};
        isC = varargin{3};

        node = [];
        ind = 1;
        indX = 1;
        x0 = 0;
        x = x0;

        if isC
            for i = 1:2
                for j = 1 : numY(2) + 1
                    y = (w-a) * (j-1) / numY(2) + a;
                    [node] = [node; ind, x, y];
                    ind = ind + 1;
                end
            end
            numX = 1;
            elem = zeros(numX * numY(2), 6);
        else
            while x0 <= b1
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

                x = x0 + dx;
                x0 = x;

                if x0 > b1 && x0 - dx < b1
                    x = b1;
                    x0 = x;
                end

                if indX == 5 && dx < 3
                    dx = dx * 2;
                    indX = 1;
                else
                    indX = indX + 1;
                end
            end
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
            end
            numX = size(node, 1)/(sum(numY) + 1) - 1;
            elem = zeros(numX * sum(numY), 6);
        end

        elemBou = cell(4, 1);
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
                if i == 1 %
                    elemBou{4} = [elemBou{4}; ind, 4];
                end

                %                 elem(ind, :) = [ind, 1, ind1, ind2, ind3, ind4];
                elem(ind, :) = [ind, 1, ind2, ind3, ind4, ind1];
                ind = ind + 1;
            end
        end

        nodeBou = cell(4, 1);
        nodeBou{1} = node(1:sum(numY)+1:end, :);
        nodeBou{2} = node(end-sum(numY):end, :);
        nodeBou{3} = node(sum(numY)+1:sum(numY)+1:end, :);
        nodeBou{4} = node(1:sum(numY)+1, :);

        % 绘制网格
        figure(1);
        hold on;
        for i = 1:size(elem, 1)
            nodes = node(elem(i, 3:end), 2:3);
            fill(nodes(:, 1), nodes(:, 2), 'w', 'EdgeColor', 'k');
        end
        if isC
            plot(node(:, 2), node(:, 3), 'ro', 'MarkerFaceColor', 'r');
        end
        title('single-edge notched beam (SE(B)) test');
        xlabel('X');
        ylabel('Y');
        axis equal;
end
end