function F = getFext(Mesh, BC, Quad, Control, Material, Pressure)
%GETFEXT External forces
%   F = GETFEXT(Mesh, BC, Quad) is a column vector of external forces
%   acting on each degree of freedom (size ndof x 1 in which ndof is the
%   number of degrees of freedom)
%
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:   Structure array with the following fields,
%           .nDOF:      Total number of DOFs
%           .ne:        Total number of elements in the mesh
%           .conn:      Array of element connectivity (size ne x nne)
%           .x:         Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%           .DOF:       Array of DOF indices (size nn x nsd)
%           .nDOFe:     Number of DOFs per element
%           .type:      the toplogical class of finite element; it is in
%                       the general form 'topology-#of nodes' ie a three
%                       node triangle is T3 a four node quadralateral is
%                       Q4 a 4 node tetrahedra is H4 a 27 node brick is
%                       B27 etc. Presently defined are L2, Q4, and Q9.
%           .nsd:       Number of spatial dimensions
%
%   BC:     Structure array with the following fields,
%           .b                          Anonymous function of distributed
%                                       body force (size 1 x nsd)
%           .traction_force_node        Column vector of nodes with
%                                       prescribed tractions
%           .traction_force_value       Column vector of prescribed tractions
%                                       on nodes
%           .traction_force_dof         Column vector of degrees of freedom
%                                       with prescribed tractions
%           .traction_force_dof_value   Column vector of prescribed tractions
%                                       on DOF
%
%   Quad:   Structure array with the following fields,
%           .W:         Vector of quadrature weights (size nq x 1)
%           .Q:         Vector of quadrature points (size nq x nsd)
%           .nq:        Number of quadrature points
%           .Nq:        Cell array (size nq x 1) with shape functions
%                       evaluated at each quadrature point
%           .dNdxiq:    Cell array (size nq x 1) with derivative of shape
%                       functions w.r.t. parent coordinates evaluated at
%                       each quadrature point

% Acknowledgements: Chris Ladubec

%% Initialize values

% initialize source (body) force vector
F = zeros(Mesh.nDOF, 1);

%% loop through all elements

for e = 1:Mesh.ne

    %% Element variables
    % nodal ids of the element's nodes
    enodes = Mesh.conn(e,:);
    % global coordinates of the element's nodes
    nodesI = Mesh.nodes(enodes,:);
    % DOFs of element nodes
    dofE = Mesh.DOF(enodes,:);
    dofE = reshape(dofE',Mesh.nDOFe,[]);
    % number of degrees of freedom
    ndofE = length(dofE);
    % elemental material
    nMat = Mesh.eMat(e,1);
    % max and min global coordinate of element
    zmax = max(nodesI(:,2));
    zmin = min(nodesI(:,2));

    if (zmax <= Mesh.Hmax) && (zmin >= Mesh.Hmin)

        %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        Q = Quad.Q;
        nq = Quad.nq;

        %% Calculate element body force

        % initialize element body force vector
        Fbe = zeros(ndofE, 1);

        % length of element. used to check that quadrature points and weights
        % are correct. A = Sum Wi*Ji
        A = 0;

        for q = 1:nq

            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};
            Nv = getNv(N, Mesh.nDOFn);

            % quadrature point in physical coordinates
            nodesi = nodesI'*N;
            r = nodesi(1,1);
            z = nodesi(2,1);

            % Jacobian of the transformation between parent and global
            % coordinates
            Je = dNdxi'*nodesI;

            % determinant of the Jacobian
            dJe = det(Je);

            % body force value at Gauss point
            % if Control.Manufactured == 1
            %     be = BC.b(r,z);
            % else
            %     be = Material.Prop(nMat).b;
            % end

            be = Material.Prop(nMat).b;

            % Applied body force
            Fbe = Fbe + 2*pi*W(q)*Nv*be*dJe*nodesi(1,1);

            % quadrature debug tool
            A = A + W(q)*dJe;
        end


        %% Calculate elemental force contribution from pore pressure

        % initialize element pore pressure force vector
        Fpe = zeros(ndofE, 1);

        % length of element. used to check that quadrature points and weights
        % are correct. A = Sum Wi*Ji
        A = 0;


        for q = 1:nq

            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};

            % quadrature point in physical coordinates
            nodesi = nodesI'*N;

            % Jacobian of the transformation between parent and global
            % coordinates
            Je = dNdxi'*nodesI;

            % determinant of the Jacobian
            dJe = det(Je);
            if dJe < 0
                error('Element %d has a negative Jacobian.', e)
            end

            % derivative of shape function in physical coordinates
            % (tensor form)
            dNdxi = dNdxi';
            B = Je\dNdxi;

            % convert B matrix to Voigt form
            r = nodesi(1,1);
            z = nodesi(2,1);
            Bv = getBv(B,N,r,Mesh.type);

            % Mapping Vector
            m = [1;1;1;0];

            % pressure at Gauss point
            % if Control.Manufactured == 1
            %     pe = Pressure.p(r,z);
            % else
                pe = Pressure.p(r) * (z <= Mesh.Hmax) * (z >= Mesh.Hmin);

                % if pe > 0
                %     pause(1)
                % end
            %end

            % Calculate local stiffness matrix
            Fpe = Fpe + 2*pi*W(q)*Bv'*m*pe*dJe*r;

            % quadrature debug tool
            A = A + W(q)*dJe;
        end



        %% Calculate element traction force vector
        % initialize element traction force vector
        Fte = zeros(ndofE, 1);
        %
        %----  Saving code for possible future adaption to edge element
        %integration-----
        %
        % 	    % loop through all tractions
        %         for i = 1:length(traction_force_node)
        %
        % 	        % check whether traction is applied to the element
        % 	        if ~isempty(find(enodes == BC.traction_force_node(i), 1)) ...
        % 	                && t_count(i) == 0
        %
        % 	            % local node number at traction location
        % 	            t_node = find(enodes == BC.traction_force_node(i),1);
        % 	            % parent coordinates of traction location
        % 	            xi = getXI(t_node,Mesh.type);
        % 	            % Shape function at traction location
        %                 [N,dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
        %                 Nv = getNv(N, Mesh.nsd);
        %
        % 	            Fte = Fte + Nv*traction_value(i,:)';
        %
        % 	            t_count(i) = 1;
        % 	        end
        %     	end

        %% Assemble element forces
        F(dofE) = F(dofE) + Fbe + Fte + Fpe;

    end
end

%% Add forces prescribed at dofs
F(BC.traction_force_dof) = F(BC.traction_force_dof) ...
    + BC.traction_force_dof_value;

end
