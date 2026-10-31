function dt = ComputeTimeStep(Mesh, Quad, Material, pbot, Control, BC)

dr = Mesh.L/Mesh.ner;
a = zeros(Mesh.ne,1);

% for each element, compute element advection matrix and add to global
for e = 1:Mesh.ne

    %% Element variables
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);    
        % global coordinates of the element's nodes
        rI = Mesh.r(enodes,:);  
        % DOFs of element nodes
        dofE = Mesh.DOF(enodes,:);
        dofE = reshape(dofE',Mesh.nDOFe,[]);
        
        % Pressure at Bottom
        pbote = pbot(dofE);

   %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        nq = Quad.nq;

        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % loop through all quadrature points
        for q = 1:nq     

            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};

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

            % Flow Properties 
            kx = Material.kx;
            muB = Material.muB;
            Sb = Material.sB;
            phi = Material.phi;

            % Calculate Advection Velocity
            dp0G = B*pbote;
            ae = (-kx/(phi*(1-Sb)*muB))*(dp0G);
            a(e) = max(a(e),abs(ae));

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end

end
    
dt = Control.Courant*dr/max(abs(a));

end



