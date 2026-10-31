function Fp2 = getFp2(Mesh, Quad, BC, t)
%   GETFP2 -> INJECTION FLUX VECTOR
%   Fp2 = getFp2(Mesh, Quad, hc, Material) is a column vector of injection 
%   flux vector for a axisymmetric vertically averaged multi-phase flow model.
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

% Acknowledgements: Giovanni Zingaro
  
%% Initialize Buoyancy Flux Vector
Fp2 = zeros(Mesh.nDOF, 1);

%% Loop through Elements

for e = 1:Mesh.ne

    % Element variables
	    % nodal ids of the element's nodes
	    enodes = Mesh.conn(e,:);    
	    % global coordinates of the element's nodes
	    rI = Mesh.r(enodes,:);  
	    % DOFs of element nodes
        dofE = Mesh.DOF(enodes,:);
        dofE = reshape(dofE',Mesh.nDOFe,[]);
        % number of degrees of freedom
        ndofE = length(dofE);

    % Shape functions and derivatives in parent coordinates
        W = Quad.W;
        nq = Quad.nq;

    % Calculate Injection Flux Vector

	    % Initialize Elemental Injection Flux Vector
	    Fp2e = zeros(ndofE, 1);

        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;
   
        % Quadature Loop
        for q = 1:nq

            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdri = Quad.dNdxiq{q};
            Nv = getNv(N, Mesh.nDOFn);

            % quadrature point in physical coordinates
            Ri = rI'*N;

            % Jacobian of the transformation between parent and global
            % coordinates
            Je = dNdri'*rI;

            % determinant of the Jacobian
            dJe = det(Je);

            % Applied body force
            Fp2e = Fp2e + W(q)*Nv*BC.s(Ri,t)*Ri*dJe;

            % quadrature debug tool
            A = A + W(q)*dJe;
        end

        Fp2(dofE,1) = Fp2(dofE,1) + Fp2e; 

end
    
end
