function Fp3 = getFp3(Mesh, Quad, hc, Material)
%   GETFP3 -> BUOYANCY FLUX VECTOR
%   Fp3 = getFp3(Mesh, Quad, hc, Material) is a column vector of buoyancy 
%   flux vector for a axisymmetric vertically averaged multi-phase flow model.
%   
% Acknowledgements: Giovanni Zingaro
  
%% Initialize Buoyancy Flux Vector and Properties
Fp3 = zeros(Mesh.nDOF, 1);

% Flow Properties
kx = Material.kx;                       % Permeability [m^2]
muC = Material.muC;                     % Viscosity of CO2 [Ns/m^2]
drho = Material.rhoB - Material.rhoC;   % Δρ = ρB - ρCO2 [kg/m^3]
g = Material.g;                         % Gravitational Constant [m/2^2]
H = Mesh.H;

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

    % Calculate Buoyancy Flux Vector

	    % Initialize Elemental Buoyancy Flux Vector
	    Fp3e = zeros(ndofE, 1);

        % Length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % Elemental Saturation of CO2
        hce = hc(dofE,1);                       % Elemental Saturation Values
   
        % Quadature Loop
        for q = 1:nq

                % Shape Functions and Derivatives in Parent Coordinates
                N = Quad.Nq{q};
                dNdxi = Quad.dNdxiq{q};
                
                % Quadrature Point in Physical Coordinates
                Ri = rI'*N;

                % Jacobian of the Transformation between Parent and Global 
                % Coordinates
                Je = dNdxi'*rI;
               
                % Determinant of the Jacobian
                dJe = det(Je);
                if dJe < 0
                    error('Element %d has a negative Jacobian.', e)
                end

                % Derivative of Shape Function in Physical Coordinates
                % (tensor form)
                dNdxi = dNdxi';
                B = Je\dNdxi;

                % Convert B Matrix to Voigt f=Form
                Bv = getBv(B', Mesh.nsd);

                % Derivative of Saturation and Saturation of CO2
                hcG = N'*hce;
                dhcG = B*hce;

                % Relative Permeability at Gauss Points
                krelC = hcG/H;

                % Elemental Force Parameters
                D = drho*hcG*krelC*kx*g/muC*dhcG;
              
                % Elemental Buoyancy Flux Vector 
                Fp3e = Fp3e + W(q)*Bv*D*Ri*dJe;

                % Quadrature Debug Tool
                A = A + W(q)*dJe; 
        end

        Fp3(dofE,1) = Fp3(dofE,1) + Fp3e; 

end
    
end
