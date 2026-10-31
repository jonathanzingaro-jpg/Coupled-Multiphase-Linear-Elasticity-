function [m_L2,m_e, d_exact, e_exact, s_exact] = ManufacturedSolutionCheck(d,s,e,Mesh, Control, Material)
%MANUFACTUREDSOLUTION_CHECK Calculates the convergence rates
%   [m_L2, m_e] = ManufacturedSolution_check(d1, d2, d3, s1, s2, s3, Mesh1,
%   Mesh2, Mesh3) calculates the rates of convergence of the L2 error norm
%   and the energy norm for using a manufactured solution
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d1, d2, d3:             Displacement vectors from 3 runs with increasing mesh
%                           refinement
%   s1, s2, s3:             Nodal stress data from 3 runs with increasing mesh
%                           refinement
%   e1, e2, e3:             Nodal strain data from 3 runs with increasing mesh
%                           refinement
%   Mesh1, Mesh2, Mesh3:    Mesh data structure from 3 runs with increasing mesh
%                           refinement
%   1 - coarsest mesh, 2 - medium mesh, 3 - finest mesh

% Manufactured solution
%   ux = x^5 + x*y^3 - y^6
%   uy = x^5 + x*y^3 - y^6
%   exx = 5x^4 + y^3
%   eyy = 3xy^2 - 6y^5
%   exy = (3xy^2 - 6y^5 + 5x^4 + y^3) 
%   sxx = E/(1+v)/(1-2v)  * [(1-v)*exx + v*eyy ]
%   syy = E/(1+v)/(1-2v)  * [ v*exx    + (1-v)*eyy ]
%   sxy = E/(1+v)/(1-2v)  * (1-2v)/2 * exy

% Acknowledgements: Jonathan Zingaro

    % Exact Solution
    r = Mesh.nodes(:,1);
    z = Mesh.nodes(:,2);

    % Exact Solution for Displacement Field
    d_exact = zeros(2*Mesh.nn,1);
    d_exact(1:2:end) = r.^6 - z.^6;
    d_exact(2:2:end) = r.^6 - z.^6;

    % Exact Solution for Stress and Strain Field
    e_exact = zeros(4, Mesh.nn);
    s_exact = zeros(4, Mesh.nn);

    e_exact(1,:) = 6.*r.^5;
    e_exact(2,:) = -6.*z.^5;
    e_exact(3,:) = (r.^6 - z.^6)./r;
    e_exact(4,:) = 6.*r.^5 - 6.*z.^5;

    E = Material.Prop(1).E; % Young's Modulus
    nu = Material.Prop(1).nu; % Poisson's Ratio

    s_exact(1,:) = -(2.*(- 19.*r.^6 + 6.*r.*z.^5 + z.^6))./(5.*r);
    s_exact(2,:) = -(2.*(- 7.*r.^6 + 18.*r.*z.^5 + z.^6))./(5.*r);
    s_exact(3,:) = -(6.*(- 3.*r.^6 + 2.*r.*z.^5 + z.^6))./(5.*r);
    s_exact(4,:) = (12.*r.^5)./5 - (12.*z.^5)./5;


    % Calculate error norms
    Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);
    nq = Quad.nq;
    
    eL2_num = 0;
    eL2_den = 0;
    eEN_num = 0;
    eEN_den = 0;
    
    for i = 1:Mesh.ne
        enodes = Mesh.conn(i,:);
        nodesI = Mesh.nodes(enodes,:);
        
        for p = 1:nq
            N = Quad.Nq{p}';
            dNdxi = Quad.dNdxiq{p};
            wp = Quad.W(p);
            Je = dNdxi'*nodesI;

            % quadrature point in physical coordinates
            nodesi = nodesI'*N';
            r = nodesi(1,1);
            
            % approximated displacement at quadrature point
            uxh_p = N*d(Mesh.xdofs(enodes));
            uyh_p = N*d(Mesh.ydofs(enodes));

            % exact displacement at quadrature point
            uxe_p = N*d_exact(Mesh.xdofs(enodes));
            uye_p = N*d_exact(Mesh.ydofs(enodes));

            % L2 norm
            eL2_num = eL2_num + 2*pi*[uxh_p - uxe_p , uyh_p - uye_p] * [uxh_p - uxe_p ; uyh_p - uye_p] * wp * r * det(Je);
            eL2_den = eL2_den + 2*pi*[uxe_p , uye_p] * [uxe_p ; uye_p] * wp * r * det(Je);
            % approximated strain at quadrature point
            ehp = e(:,enodes)*N';

            % exact strain at quadrature point
            eep = e_exact(:,enodes)*N';
            D = getD(Material,Mesh, 1);

            % approximated stress at quadrature points
            shp = s(:,enodes)*N';
            sep = s_exact(:,enodes)*N';

            % e norm
            eEN_num = eEN_num + (ehp-eep)'*(shp - sep) * wp * det(Je);
            eEN_den = eEN_den + eep'*sep * wp * det(Je);
            
        end
    end

    m_L2 = sqrt(eL2_num/eL2_den);
    m_e = (eEN_num.^0.5)./(eEN_den.^0.5);
    
end


