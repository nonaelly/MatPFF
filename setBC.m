function BoundaryCondition = setBC(fixNode, nodeForce, NDof)
% SENT - Single Edge Notched Tension
% % ** code by P.M.Hu @bit.edu.cn (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com

MoveNode = fixNode(fixNode(:,3)>0, 1:2);

% Dirichlet
BoundaryCondition.DirchletDOF = fixNode(:,1)*2 + (fixNode(:,2)-2);
BoundaryCondition.Dirichlet   = fixNode(:,3);
BoundaryCondition.FreeDOF     = setdiff([1:NDof]', BoundaryCondition.DirchletDOF);
BoundaryCondition.BDforce     = MoveNode(:,1)*2 + (MoveNode(:,2)-2);

% Neumann(Natural) Boundary Condition
F = zeros(NDof, 1);
for i = 1 : size(nodeForce, 1)
    F(nodeForce(i, 1)*2-2 + nodeForce(i, 2)) = F(nodeForce(i, 1)*2-2 + nodeForce(i, 2)) + nodeForce(i, 3);
end
BoundaryCondition.RHS         = F;

end