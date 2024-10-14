function fixNode = generateBC(BC, bou)

fixNode = [];
for i = 1 : 4
    line = bou{i};
    switch BC{i}
        case {'dy', 'y'}
            fixNode = [fixNode; line(:, 1), 2*ones(size(line, 1), 1), zeros(size(line, 1), 1)];
        case {'dx', 'x'}
            fixNode = [fixNode; line(:, 1), 1*ones(size(line, 1), 1), zeros(size(line, 1), 1)];
        case {'p', 'f'}
        case {'d'}
            for j = 1 : size(line, 1)
                fixNode = [fixNode; line(j, 1), 1, 0; line(j, 1), 2, 0];
            end
    end
end
fixNode = sortrows(fixNode, 1);
end