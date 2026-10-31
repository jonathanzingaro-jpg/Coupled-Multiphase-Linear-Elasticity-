function Fp4 = getFp4(Mesh, Quad, h, Material)
%   GETFP3 -> ELEVATION FLUX VECTOR
%   Fp3 = getFp3(Mesh, Quad, hc, Material) is a column vector of elevation 
%   flux vector for a axisymmetric vertically averaged multi-phase flow model.
%   
% Acknowledgements: Giovanni Zingaro
  
%% Initialize Buoyancy Flux Vector
Fp4 = zeros(Mesh.nDOF, 1);

% Retrieve Saturation of Brine and CO2
hb = h.hb;
hc = h.hc;

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

    % Calculate Elevation Flux Vector

	    % Initialize Elemental Elevation Flux Vector
	    Fp4e = zeros(ndofE, 1);

        % Length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % Elemental Saturation of CO2 and Brine
        hce = hc(dofE,1);                       % Elemental Saturation Values
        hbe = hb(dofE,1);                       % Elemental Saturation Values

        % Flow Properties 
        kx = Material.kx;                       % Permeability [m^2]
        muC = Material.muC;                     % Viscosity of CO2 [Ns/m^2]
        rhoB = Material.rhoB;                   % Density of Brine [kg/m^3]
        g = Material.g;                         % Gravitational Constant [m/2^2]
    
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

                % Convert B Matrix to Voigt Form
                Bv = getBv(B', Mesh.nsd);

                % Derivative of Saturation and Saturation of CO2
                hcG = N'*hce;
                He = hce + hbe;
                dHe = B*He;

                % Constutitive Law
                D = hcG * kx * krelC / muC * rhoB * g * dHe;
              
                % Elemental Buoyancy Flux Vector 
                Fp4e = Fp4e + W(q)*Bv*D*Ri*dJe;

                % Quadrature Debug Tool
                A = A + W(q)*dJe; 
        end

        Fp4(dofE,1) = Fp4(dofE,1) + Fp4e; 

end
    
end
