function Kg = getKsu(Mesh, Quad, Material, p0, Control)

% getKp Diffusion Matrix for Pressure Equation used in Vertically Averaged 
%   Multi-Phased Fluid Flow Kp = getKp(Mesh, Quad, Material) returns the 
%   diffusion matrix for Vertically Averaged Multi-Phased Fluid Flow

% Created by: Jonathan Zingaro

% Initialize Stiffness Matrix 
vec_size = Mesh.ne*(Mesh.nne * Mesh.nsd)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Kvec = zeros(vec_size, 1);                  % vectorized stiffness matrix
count = 1;                                  % DOF counter

% Flow Properties
kx = Material.kx;
muB = Material.muB;
sB = Material.sB;
phi = Material.phi;

% Temporal and Spatial Discretization Properties
dr = Mesh.dr;
theta = Control.theta;
dt = Control.TimeStep;

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

        % Elemental Pressure at Bottom of Reservoir, Slope at Bottom of
        % Reservoir 
        p0e = p0(dofE,1);

        % loop through all quadrature points
        for q = 1:nq     

            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};

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
        
            % Calculate the First and Second Derivatives of Pbot
            dp0G = Bv'*p0e; % First Derivative of Pbot

            % Calculate Advection Velocity 
            ae = -(kx/(phi*(1-sB)*muB))*(dp0G);

            % Stabilization Paramater
            te = Control.Stab_Scale*(dr/(2*abs(ae)));
            
            % Calculate local stiffness matrix
            Ke = Ke - W(q)*Bv*ae*N'*Ri*dJe + W(q)*Bv*(ae^2)*te*Bv'*Ri*dJe;
            
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
Kg = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);

end



