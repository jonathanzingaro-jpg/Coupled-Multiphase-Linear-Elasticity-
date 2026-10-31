function [stress, strain] = getStrain(d, Mesh, Material, calc_type, Quad, ~)
%GETSTRAIN Evaluate stress and strain
%   [strain, stress] = GETSTRAIN(d, Mesh, Material) returns two matrices of 
%   strains and stresses computed at the center of each element. The 
%   matrices are of size dim x ne, in which dim = 1 for 1D elements, 3 for 
%   2D elements, and 6 for 3D elements.
%
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'none') does not
%   compute the stresses or strains
% 
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'nodal') returns 
%   matrices of nodal-averaged strains (size dim x nn).
% 
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'center') returns 
%   matrices of strains computed at the center of each element 
%   (size dim x ne).
%
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'L2projection') returns 
%   matrices of L2-projected stresses at the nodes (size dim x nn).
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   d:          Vector of displacements at each DOF (size ndof x 1 in 
%               which ndof is the number of degrees of freedom)
% 
%   Mesh:       Structure array with the following fields,
%               .nsd:   Number of spatial dimensions
%               .ne:    Total number of elements in the mesh
%               .nn:    Total number of nodes 
%               .nne:   Vector of number of nodes per element (size nsd x 1)
%               .type:  the toplogical class of finite element; it is in 
%                       the general form 'topology-#of nodes' ie a three 
%                       node triangle is T3 a four node quadralateral is 
%                       Q4 a 4 node tetrahedra is H4 a 27 node brick is 
%                       B27 etc. Presently defined are L2, Q4, and Q9. 
%               .nDOFe: Number of DOFs per element
%               .conn:  Array of element connectivity (size ne x nne)
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .DOF:   Array of DOF indices (size nn x nsd)
% 
%   Material:   Structure array with the following fields,
%               .E:     Modulus of elasticity
%               .nu:    Poisson's ratio
%               .Dtype: 2D approximation ('PlaneStrain' or 'PlainStress')
%
% 	Quad: 	Structure array with following fields,
%       	.nq:     Number of quadrature points	
%       	.W:      Vector of quadrature weights (size nq x 1)      
%       	.Q:      Vector of quadrature points	(size nq x nsd)
%       	.Nq:     Cell array (size nq x 1) with shape functions  
%       	         evaluated at each quadrature point
%       	.dNdxiq: Cell array (size nq x 1) with derivative of shape 
%       	         functions w.r.t. parent coordinates evaluated at 
%       	         each quadrature point
%       	.Nv:     Cell array (size nq x 1) with shape functions 
%       	         evaluated at each quadrature point in Voigt form
%
%   Acknowledgements: Jonathan Zingaro

    % Initialize Stress and Strain Arrays
    switch calc_type
        case 'Nodal'
            strain = zeros(4, Mesh.nn);
            stress = zeros(4, Mesh.nn);
            count = zeros(4, Mesh.nn);
        case 'Center'
            strain = zeros(4, Mesh.ne);
            stress = zeros(4, Mesh.ne);
    end

    % Initialize Centriod of Elements 
    switch Mesh.type
        case 'T3'
            nodesi = [1/3,1/3];
        case 'Q4'
            nodesi = [0,0];
    end

    % Loop through all elements
    for e = 1:Mesh.ne

        %% Element variables
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);
        % global coordinates of the element's nodes
        nodesI = Mesh.nodes(enodes,:);
        % DOFs of element nodes
        dofE = Mesh.DOF(enodes,:);
        dofE = reshape(dofE',Mesh.nDOFe,[]);
        % local displacement vector
        de = d(dofE);
        % elemental material
        nMat = Mesh.eMat(e,1);
        
        %% Constitutive matrix
        D = getD(Material, Mesh, nMat);

        switch calc_type
            case 'Center'

                [N, dNdxi] = lagrange_basis(Mesh.type, nodesi, Mesh.nsd);
                Je = dNdxi'*nodesI;
                B = Je\(dNdxi');
                r = mean(nodesI(:,1));
                Bv = getBv(B,N,r,Mesh.type);
                strain(:, e) = Bv*de;
                stress(:, e) = D*strain(:, e);

            case 'Nodal'
                % initialize strain element
                strain_e = zeros(4, Mesh.nne);
                stress_e = zeros(4, Mesh.nne);

                % loop through all nodes and calculate nodal strains/stresses
                for n = 1:Mesh.nne

                    % node point in parent coordinate
                    nodesi = getXI(n, Mesh.type);

                    % Shape function and derivatives in parent coordinates
                    [N, dNdxi] = lagrange_basis(Mesh.type, nodesi, Mesh.nsd);

                    % Jacobian of the transformation between parent and global
                    % coordinates
                    Je = dNdxi'*nodesI;

                    % derivative of shape function in physical coordinates
                    % (tensor form)
                    dNdxi = dNdxi';
                    B = Je\dNdxi;

                    % radial coordinate 
                    r = nodesI(n,1);

                    % convert B matrix to Voigt form
                    Bv = getBv(B,N,r,Mesh.type);

                    % strain_e = [strainx_n1  strainx_n2...;
                    %strainy_n1  strainy_n2...;
                    %strainxy_n1 strainxy_n2...];
                    strain_e(:,n) = Bv*de;
                    stress_e(:,n) = D*strain_e(:,n);
                end
                % Add to global strains
                strain(:,enodes) = strain(:,enodes) + strain_e;
                stress(:,enodes) = stress(:,enodes) + stress_e;
                count(:,enodes) = count(:,enodes) + 1;
          
        end
    end


    switch calc_type
        case 'Nodal'
            % For nodal strains, divide by count to get the average
            strain = strain./count;
            stress = stress./count;
    end
end





