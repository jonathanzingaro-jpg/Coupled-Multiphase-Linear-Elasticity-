function N = getxi(Mesh, e, xp)
%GETXI Parent coordinates of the nodal position specified.
%   xi_node = GETXI(node_num, type) returns the parent coordinates at the  
%   node position (node_num) for the specified element type (type). 
%   For multi-dimensional elements, the coordinates are given as a 
%   vector of size 1 x nsd where nsd is the number of spatial dimensions.
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   node_num:   node number, based on standard element numbering
%   type:       the topological class of finite element; it is in the 
%               general form 'topology-#of nodes' ie a three node 
%               triangle is T3 a four node quadralateral is Q4 a 4 node 
%               tetrahedra is H4 a 27 node brick is B27 etc. Presently 
%               defined are L2, L3, L4, T3, Q4, Q9, and B8 B27.  

switch Mesh.type
    case 'L2'
        
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);

        % global coordinates of the element's nodes
        xI = Mesh.r(enodes,:);  

        % Shape Function in Global Coordinates
        x1 = min(xI);
        x2 = max(xI);

        % Shape Function in Global Coordinates
        N = (1/(x2-x1))*[x2-xp,xp-x1];

    case 'L3'
        
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);

        % global coordinates of the element's nodes
        xI = Mesh.x(enodes,:);

        % Ensure that the nodes are in the correct order
        % xI should be a 3x1 vector for a three-node element
        x1 = xI(1);
        x2 = xI(2);
        x3 = xI(3);

        % Shape Function in Global Coordinates
        % Evaluate shape functions at point xp
        % Make sure xp is defined in your code; it is the position where you evaluate N
        N1 = ((xp - x2) * (xp - x3)) / ((x1 - x2) * (x1 - x3));
        N2 = ((xp - x1) * (xp - x3)) / ((x2 - x1) * (x2 - x3));
        N3 = ((xp - x1) * (xp - x2)) / ((x3 - x1) * (x3 - x2));

        % Assemble shape function vector for the element
        N = [N1, N2, N3];


end

