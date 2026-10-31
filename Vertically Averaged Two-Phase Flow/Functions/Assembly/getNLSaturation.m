function [H_hat, R_hat, Fint] = getNLSaturation(Mesh, Quad, Material, pbot, h, Control, C)
%   GETNLSATURATION computes the Hessian matrix (H) and the residual vector 
%   for a nonlinear axisymmetric vertically averaged multi-phase flow
%   problem, where the nonlinearity arises as through a linear relative
%   permeability relation for the brine phase. 
%
% Acknowledgements: Giovanni Zingaro
% -------------------------------------------------------------------------

% Initialize Hessian Matrix
vec_size = Mesh.ne*(Mesh.nne)^2;            % vector size (dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Hvec = zeros(vec_size, 1);                  % vectorized Hessian matrix
Kvec = zeros(vec_size, 1);                  % vectorized Advection matrix


count = 1;                                  % DOF counter

% Initialize Internal Flux Vector
Fint = zeros(Mesh.nDOF, 1);  

% Initialize hb Vector
hb = h.hb;      % Get hb from h structure for timestep n
hbn1 = h.hbn1;  % Get hbn1 from h structure for timestep n-1

% Flow Properties for Brine Phase
kr = Material.kx;       % Radial Permeablity of Brine
muB = Material.muB;     % Viscoisty of Brine
sB = Material.sB;       % Residual Saturation of Brine
phi = Material.phi;     % Porosity of Rock Mass
H = Mesh.H;             % Height of Reservoir

% Temporal and Spatial Discretization Properties
dr = Mesh.dr;
theta = Control.theta;
dt = Control.TimeStep;

% For each element, compute element hessian matrix and add to global
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
        % saturation at element's nodes
        hbe = hb(dofE);

    %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        nq = Quad.nq;

    %% Assemble hessian matrix
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element hessian matrix
        He = zeros(ndofE, ndofE); 

        % initialize element advection matrix
        Ke = zeros(ndofE, ndofE); 
        
        % initialize element internal flux vector
        fe = zeros(ndofE, 1); 

        % Elemental Pressure at Bottom of Reservoir, Slope at Bottom of
        % Reservoir
        pbote = pbot(dofE,1);

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

            % Compute Saturation at Gaussian Point
            hbG = N'*hbe;

            % Calculate the dPbot/dr
            dp0G = Bv'*pbote; 

            % Calculate Advection Velocity at Gauss Point
            krelG = hbG/H;
            ae = -(kr*krelG/(phi*(1-sB)*muB))*(dp0G);
            dadr = -(kr/(phi*(1-sB)*muB*H))*(dp0G);

            % Calculate Elemental Hessian Matrix
            He = He ... 
                 - W(q)*Bv*ae*N'*Ri*dJe;% ...
                 %- W(q)*Bv*dadr*N'*hbe*N'*Ri*dJe;

            % Assmeble Elemental Advection Matrix
            Ke = Ke - W(q)*Bv*ae*N'*Ri*dJe;
            
            % Calculate local internal flux vector
            fe = fe - W(q)*Bv*ae*N'*Ri*dJe*hbe;

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized stiffness matrix
        count = count + ndofE^2;
        He = reshape(He, [ndofE^2, 1]);
        Ke = reshape(Ke, [ndofE^2, 1]);
        rowmatrix = dofE*ones(1,ndofE);        
        rowe = reshape(rowmatrix,[ndofE^2, 1]);
        cole = reshape(rowmatrix',[ndofE^2, 1]);
        
        Hvec(count-ndofE^2:count-1) = He;
        Kvec(count-ndofE^2:count-1) = Ke;
        row(count-ndofE^2:count-1) = rowe;
        col(count-ndofE^2:count-1) = cole;
       
        %% Assemble element internal forces
        Fint(dofE) = Fint(dofE) + fe;
end

% Sparse Hessian Matrix
H = sparse(row, col, Hvec, Mesh.nDOF, Mesh.nDOF);

% Sparse Advection Matrix
K = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);

% Hessian matrix in transient case
H_hat = (theta*H + C./dt);

% Residual
R_hat = (C + dt*theta*K) * hb - (C)*hbn1;

end

