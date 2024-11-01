% -------------------------------- %
% Copyright (c) 2021, Carlos Souto %
% csouto@fe.up.pt                  %
% -------------------------------- %

function parts = loadinp(inp)
    
    % ------------------------------------------ %
    % Loads the inp data to a parts struct array %
    % ------------------------------------------ %
    
    % Get and process text from the file
    text = fileread(inp);         % read all the text
    text = strtrim(text);         % remove leading and trailing whitespaces
    text = strrep(text, ' ', ''); % remove any in-between whitespaces
    text = lower(text);           % lower the case
    text = inpremcom(text);       % remove comment lines from inp text data
    
    % Split text data in parts
    data = inpsplit(text);
    numparts = size(data, 1);
    
    % Extract data
    for i = 1:1:numparts
        parts(i).Name = getname(data{i});
        parts(i).Nodes = getnodes(data{i});
        parts(i).Elements = getelements(data{i});
    end
    
end

function output = inpremcom(text) 
    
    % ---------------------------------------- %
    % Removes comment lines from inp text data %
    % ---------------------------------------- %
    
    text = splitlines(text);
    
    output = {};
    j = 1;
    
    for i = 1:1:size(text, 1)
        if ~startsWith(text{i}, '**')
            output{j, 1} = text{i};
            j = j + 1;
        end
    end
    
    output = strjoin(output, '\n');
    
end

function output = inpsplit(text)
    
    % ----------------------------------- %
    % Splits the inp text data into parts %
    % ----------------------------------- %
    
    a = strsplit(text, '*part')';
    b = cell(size(a, 1) - 1, 1);
    for i = 1:1:size(b, 1)
        c = strsplit(a{i + 1}, '*endpart');
        b{i} = strtrim(c{1});
    end
    output = b;
    
end

function output = getname(partdata)
    
    % -------------------------------- %
    % Gets the name from the part data %
    % -------------------------------- %
    
    lines = splitlines(partdata);
    firstline = lines{1};
    split = strsplit(firstline, '=');
    output = split{2};
    
end

function output = getnodes(partdata)
    
    % --------------------------------- %
    % Gets the nodes from the part data %
    % --------------------------------- %
    
    lines = splitlines(partdata);
    for i = 1:1:size(lines, 1)
        line = lines{i};
        if contains(line, '*node')
            for j = (i + 1):1:size(lines, 1)
                line = lines{j};
                if contains(line, '*')
                    return
                else
                    data = str2num(line);
                    output(data(1), :) = data(2:end);
                end
            end
        end
    end
    
end

function output = getelements(partdata)
    
    % ------------------------------------ %
    % Gets the elements from the part data %
    % ------------------------------------ %
    
    lines = splitlines(partdata);
    for i = 1:1:size(lines, 1)
        line = lines{i};
        if contains(line, '*element')
            split = strsplit(line, '=');
            eltype = upper(split{2});
            for j = (i + 1):1:size(lines, 1)
                line = lines{j};
                if contains(line, '*')
                    break
                else
                    data = str2num(line);
                    element.Type = eltype;
                    element.Connectivity = data(2:end);
                    output(data(1), 1) = element;
                end
            end
        end
    end
end
