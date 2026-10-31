function Cg = getCsu(Mesh, Quad)
%   getCg assembles the unstablized storage matrix (Cg) for the saturation 
%   equation used in the formulation of 1D multi-phased flow model using
%   axisymmetric coordinates. This function assembles the LSFEM
%   stabilization matrix.
% ------------------------------------------------------------------------
% Created by: Jonathan Zingaro

% Initialize Storage Matrix 
vec_size = Mesh.ne*(Mesh.nne * Mesh.nsd)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Cvec = zeros(vec_size, 1);                  % vectorized stiffness matrix
count = 1;                                  % DOF counter

% for each element, compute the elemental storage matrix for
% the saturation equation for Vertically Averaged Multi-phased flow. 
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

    %% Assemble Storage Matrix
    
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element storage matrix
        Ce = zeros(ndofE, ndofE); 

        % loop through all quadrature points
        for q = 1:nq     

            % Shape functions, derivatives and location of Gauss Point in parent coordinates
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

            % Calculate local storage matrix
            Ce = Ce + W(q)*(N)*N'*Ri*dJe; 

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized storage matrix
        count = count + ndofE^2;
        Ce = reshape(Ce, [ndofE^2, 1]);
        rowmatrix = dofE*ones(1,ndofE);        
        rowe = reshape(rowmatrix,[ndofE^2, 1]);
        cole = reshape(rowmatrix',[ndofE^2, 1]);
        
        Cvec(count-ndofE^2:count-1) = Ce;
        row(count-ndofE^2:count-1) = rowe;
        col(count-ndofE^2:count-1) = cole;
end

% sparse storage matrix
Cg = sparse(row, col, Cvec, Mesh.nDOF, Mesh.nDOF);

end
