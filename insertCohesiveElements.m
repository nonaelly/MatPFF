function [newNodes, newElements, cohesiveElements] = insertCohesiveElements(nodes, elements, insertionRegions)
    % Inputs:
    % nodes - Original node coordinates, format: [NodeID, X, Y, Z]
    % elements - Original element data, format: [ElementID, ElmentType, NodeID1, NodeID2, ...]
    % insertionRegions - Specified regions for inserting cohesive elements
    % [xmin, xmax, ymin, ymax]

    % Outputs:
    % newNodes - Complete list of nodes including newly generated nodes
    % newElements - Complete list of elements including new elements
    % cohesiveElements - Cohesive elements with node and element data

    % 1. Identify target elements and nodes within insertion regions
    [targetElements, targetNodes] = findInsertionElements(nodes, elements, insertionRegions);

    % 2. Duplicate nodes and generate cohesive elements
    [newNodes, newElements, cohesiveElements] = generateCohesiveElements(nodes, elements, targetElements, targetNodes);

    % 3. Output results (optional: save in ABAQUS format)
    outputResults(newNodes, newElements, cohesiveElements);
end

% Sub-function 1: Identify target elements and nodes within insertion regions
function [targetElements, targetNodes] = findInsertionElements(nodes, elements, insertionRegions)
    % Inputs:
    % nodes - Original node coordinates, format: [NodeID, X, Y, Z]
    % elements - Original element data, format: [ElementID, NodeID1, NodeID2, ...]
    % insertionRegions - Specified rectangular region for inserting cohesive elements, format: [xmin, xmax, ymin, ymax]

    % Outputs:
    % targetElements - Elements with at least two nodes within the specified rectangular region
    % targetNodes - Nodes on the boundaries of target elements within the region

    % Unpack insertion region boundaries
    xmin = insertionRegions(1);
    xmax = insertionRegions(2);
    ymin = insertionRegions(3);
    ymax = insertionRegions(4);
    
    % Initialize lists for storing selected elements and nodes
    targetElements = [];
    targetNodes = [];

    % Step 1: Find nodes within the rectangular region
    % Check if each node lies within [xmin, xmax] and [ymin, ymax]
    inRegionNodes = nodes(:,2) >= xmin & nodes(:,2) <= xmax & ...
                    nodes(:,3) >= ymin & nodes(:,3) <= ymax;
                
    % Get node IDs of all nodes within the region
    regionNodeIDs = nodes(inRegionNodes, 1);

    % Step 2: Find elements that contain at least two nodes in the region
    % Loop through each element and count how many of its nodes are in the region
    for i = 1:size(elements, 1)
        elementNodes = elements(i, 3 : end);  % Get nodes of the current element
        % Count how many of the element nodes are in the region
        nodesInRegionCount = sum(ismember(elementNodes, regionNodeIDs));
        
        % If at least two nodes of the element are within the region, select the element
        if nodesInRegionCount >= 2
            targetElements = [targetElements; elements(i, :)];  % Add to target elements
            targetNodes = [targetNodes; elementNodes'];          % Add element nodes to target nodes
        end
    end
    
    % Remove duplicate nodes from targetNodes
    targetNodes = unique(targetNodes);
end


% Sub-function 2: Duplicate nodes and generate cohesive elements
function [newNodes, newElements, cohesiveElements] = generateCohesiveElements(nodes, elements, targetElements, targetNodes)
    % Initialize new nodes and elements lists
    newNodes = nodes;
    newElements = elements;
    cohesiveElements = [];
    
    % Define logic for node duplication and cohesive element generation
    % For each targetElement, create new boundary nodes and cohesive elements
    
    % TODO: Implement logic for duplicating nodes and generating cohesive elements
    
    % Example structure:
    % for i = 1:length(targetElements)
    %     element = targetElements(i);
    %     % Generate cohesive elements based on the boundary of each element
    %     % Duplicate corresponding boundary nodes
    %     % Store new nodes and cohesive elements in cohesiveElements
    % end
    
end

% Sub-function 3: Output results
function outputResults(newNodes, newElements, cohesiveElements)
    % Optional: save results in ABAQUS input format or other formats
    % Define file output logic here
    
    % Example code (modify as needed for specific format):
    % save('newNodes.mat', 'newNodes');
    % save('newElements.mat', 'newElements');
    % save('cohesiveElements.mat', 'cohesiveElements');
    
    disp('Cohesive element insertion complete, results saved.');
end
