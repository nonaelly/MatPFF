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

% Step 1: Find nodes within the rectangular region
% Check if each node lies within [xmin, xmax] and [ymin, ymax]
inRegionNodes = nodes(:,2) >= xmin & nodes(:,2) <= xmax & ...
    nodes(:,3) >= ymin & nodes(:,3) <= ymax;

% Get node IDs of all nodes within the region
regionNodeIDs = nodes(inRegionNodes, 1);
targetNodes = regionNodeIDs;

% Step 2: Find elements that contain at least two nodes in the region
% Loop through each element and count how many of its nodes are in the region
for i = 1:size(elements, 1)
    elementNodes = elements(i, 3 : end);  % Get nodes of the current element
    % Count how many of the element nodes are in the region
    nodesInRegionCount = sum(ismember(elementNodes, regionNodeIDs));

    % If at least two nodes of the element are within the region, select the element
    if nodesInRegionCount >= 2
        [targetElements] = [targetElements; elements(i, :)];  % Add to target elements
    end
end

end

% Sub-function 2: Duplicate nodes and generate cohesive elements
function [newNodes, newElements, cohesiveElements] = generateCohesiveElements(nodes, elements, targetElements, targetNodes)
% Inputs:
% nodes - Original node coordinates, format: [NodeID, X, Y, Z]
% elements - Original element data, format: [ElementID, ElementType, NodeID1, NodeID2, ...]
% targetElements - Elements in which cohesive elements are to be inserted
% targetNodes - Nodes that are within the specified insertion region

% Outputs:
% newNodes - Updated list of nodes, including newly generated nodes for cohesive elements
% newElements - Updated list of elements, including new cohesive elements
% cohesiveElements - Newly created cohesive elements with node and element data

% Initialize outputs and data structures
newNodes = nodes;
newElements = elements;
cohesiveElements = [];
edgeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');  % To track shared edges
% Temporary list to store the new nodes created for cohesive elements
newNodesTemp = [];

idxCo = size(elements, 1);
% Loop through target elements to identify shared edges and insert cohesive elements
for i = 1:size(targetElements, 1)
    elementID = targetElements(i, 1);
    elementNodes = targetElements(i, 3:end);  % Nodes of the current element

    % Loop over each edge of the element
    for j = 1:length(elementNodes)
        % Get the node IDs of the current edge (each edge has two nodes)
        nodeA = elementNodes(j);
        nodeB = elementNodes(mod(j, length(elementNodes)) + 1);  % Wrap around to form edges

        % Skip edges where both nodes are not within the insertion region
        if ~ismember(nodeA, targetNodes) || ~ismember(nodeB, targetNodes)
            continue;
        end

        % Sort node IDs to create a consistent edge key
        edgeKey = sprintf('%d-%d', min(nodeA, nodeB), max(nodeA, nodeB));

        if isKey(edgeMap, edgeKey)
            % Edge is already mapped, meaning it is shared between two elements
            existingData = edgeMap(edgeKey);

            % Retrieve the nodes and data from the other element that shares this edge
            newID_A = existingData.newNodes(1);
            newID_B = existingData.newNodes(2);

            % Calculate centroids of left and right elements
            leftNodesCoords = nodes(existingData.leftNodes, 2:3);
            leftCentroid = mean(leftNodesCoords, 1);

            rightElementNodes = elementNodes;  % All nodes of the right element (current element)
            rightNodesCoords = nodes(rightElementNodes, 2:3);
            rightCentroid = mean(rightNodesCoords, 1);

            % Determine if we should place nodes in counterclockwise order
            idxCo = idxCo + 1;
            if isCounterClockwise(leftCentroid, rightCentroid, nodes(nodeA, 2:3), nodes(nodeB, 2:3))
                % Counter-clockwise order
                cohesiveElements = [cohesiveElements; [idxCo, 2, nodeA, nodeB, newID_A, newID_B]];
            else
                % Adjust to counter-clockwise order
                cohesiveElements = [cohesiveElements; [idxCo, 2, newID_A, newID_B, nodeA, nodeB]];
            end

            newElements(elementID, 2 + [j, mod(j, length(elementNodes)) + 1]) = [newID_B, newID_A];
        else
            % Edge is not shared, so we duplicate nodes for the cohesive element

            if isempty(newNodesTemp)
                nodeACoord = nodes(nodeA, 2:3);  % Coordinates of node A
                [newNodesTemp] = [newNodesTemp; size(newNodes, 1) + 1, nodeACoord];  % Add node A to newNodesTemp
                newID_A = size(newNodes, 1) + 1;  % New ID for node A
                % Now add these new nodes to the global newNodes list
                [newNodes] = [newNodes; newID_A, nodeACoord];
                % Check for node B
                nodeBCoord = nodes(nodeB, 2:3);  % Coordinates of node B
                [newNodesTemp] = [newNodesTemp; size(newNodes, 1) + 1, nodeBCoord];  % Add node B to newNodesTemp
                newID_B = size(newNodes, 1) + 1;  % New ID for node B
                % Now add these new nodes to the global newNodes list
                [newNodes] = [newNodes; newID_B, nodeBCoord];
            else
                % Add new cohesive nodes only if they are not already present
                % Check for node A
                nodeACoord = nodes(nodeA, 2:3);  % Coordinates of node A
                if ~any(ismember(newNodesTemp(:, 2:3), nodeACoord, 'rows'))  % If node A is not in newNodesTemp
                    [newNodesTemp] = [newNodesTemp; size(newNodes, 1) + 1, nodeACoord];  % Add node A to newNodesTemp
                    newID_A = size(newNodes, 1) + 1;  % New ID for node A
                    % Now add these new nodes to the global newNodes list
                    [newNodes] = [newNodes; newID_A, nodeACoord];
                else
                    % Node A already exists, so find its ID from newNodes
                    newID_A = newNodesTemp(ismember(newNodesTemp(:, 2:3), nodeACoord, 'rows'));  % Find existing node A ID
                end

                % Check for node B
                nodeBCoord = nodes(nodeB, 2:3);  % Coordinates of node B
                if ~any(ismember(newNodesTemp(:, 2:3), nodeBCoord, 'rows'))  % If node B is not in newNodesTemp
                    [newNodesTemp] = [newNodesTemp; size(newNodes, 1) + 1, nodeBCoord];  % Add node B to newNodesTemp
                    newID_B = size(newNodes, 1) + 1;  % New ID for node B
                    % Now add these new nodes to the global newNodes list
                    [newNodes] = [newNodes; newID_B, nodeBCoord];
                else
                    % Node B already exists, so find its ID from newNodes
                    newID_B = newNodesTemp(ismember(newNodesTemp(:, 2:3), nodeBCoord, 'rows'));  % Find existing node B ID
                end

            end

            % Store edge data in edgeMap for future reference by shared elements
            edgeMap(edgeKey) = struct('newNodes', [newID_A, newID_B], 'leftNodes', elementNodes);
        end
    end
end
end

% Helper function to check if the order of nodes is counter-clockwise
function ccw = isCounterClockwise(leftCentroid, rightCentroid, nodeA, nodeB)
% This function returns true if the cohesive element formed by the edge
% should be arranged in a counter-clockwise order based on the centroids of the left and right elements.

vector1 = rightCentroid - leftCentroid;  % Vector from left centroid to right centroid
vector2 = nodeB - nodeA;                 % Vector along the edge from nodeA to nodeB
crossProd = vector1(1) * vector2(2) - vector1(2) * vector2(1);
ccw = crossProd > 0;  % Positive cross product indicates counter-clockwise order
end


