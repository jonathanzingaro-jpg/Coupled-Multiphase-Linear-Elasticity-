function Kp = getKp(Mesh, Quad, Material, h)
%   getKp Diffusion Matrix for Pressure Equation used in Vertically Averaged 
%   Multi-Phased Fluid Flow Kp = getKp(Mesh, Quad, Material) returns the 
%   diffusion matrix for Vertically Averaged Multi-Phased Fluid Flow in
%   axisymmetric coordinates.

% Created by: Jonathan Zingaro

% Initialize Stiffness Matrix 
vec_size = Mesh.ne*(Mesh.nne * Mesh.nsd)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Kvec = zeros(vec_size, 1);                  % vectorized stiffness matrix
count = 1;                                  % DOF counter

% Flow Properties
kx = Material.kx;                           % Permability [m^2]
muC = Material.muC;                         % Viscosity of CO2 [Ns/m^2]
muB = Material.muB;                         % Viscosity of Brine [Ns/m^2]
H = Mesh.H;                                 % Height of Reservoir

% for each element, compute the elemental diffusion/stiffness matrix for
% the pressure equation for Vertically Averaged Multi-phased flow. 
for e = 1:Mesh.ne

    %% Element variables
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);    
        % global coordinates of the element's nodes
        rI = Mesh.r(enodes,:);  
        % DOFs of element nodes
        dofE = Mesh.DOF(enodes,:);
        dofE = reshape(dofE',Mesh.nDOFe,[]);
        % number of degrees of freedom
        ndofE = numel(dofE);
               
    %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        nq = Quad.nq;

    %% Assemble Diffusion/Stiffness Matrix
    
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element stiffness matrix
        Ke = zeros(ndofE, ndofE); 

        % Elemental Saturations for both CO2 and Brine Phase
        hce = h.hcn1(dofE);
        hbe =  h.hbn1(dofE);

        % loop through all quadrature points
        for q = 1:nq     

            % Shape functions, derivatives and location of Gauss Point in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};
            rG = Quad.Q(q); 

            % quadrature point in physical coordinates
            Ri = rI'*N;

            % Jacobian of the transformation between parent and global 
            % coordinates
            Je = dNdxi'*rI;
           
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
            Bv = getBv(B', Mesh.nsd);

            % Saturations at Quadrature Points
            hbG = N'*hbe;
            hcG = N'*hce;

            % Relative Permeability at Gauss Points
            krelB = hbG/H;
            krelC = hcG/H;

            % Constutitive Law 
            D = ((hcG*krelC/muC)+(hbG*krelB/muB))*kx;

            % Calculate local stiffness matrix
            Ke = Ke + W(q)*Bv*D*Bv'*Ri*dJe;

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized stiffness matrix
        count = count + ndofE^2;
        Ke = reshape(Ke, [ndofE^2, 1]);
        rowmatrix = dofE*ones(1,ndofE);        
        rowe = reshape(rowmatrix,[ndofE^2, 1]);
        cole = reshape(rowmatrix',[ndofE^2, 1]);
        
        Kvec(count-ndofE^2:count-1) = Ke;
        row(count-ndofE^2:count-1) = rowe;
        col(count-ndofE^2:count-1) = cole;
end

% sparse stiffness matrix
Kp = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);

end
