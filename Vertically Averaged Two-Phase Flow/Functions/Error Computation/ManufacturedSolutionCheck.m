function [m_L2, hb_exact] = ManufacturedSolutionCheck(hb,Mesh, Control, t)
%   MANUFACTUREDSOLUTION_CHECK Calculates the convergence rates
%   [m_L2, m_e] = ManufacturedSolution_check(hb,Mesh, Control, t)
%   calculates the rates of convergence of the L2 error norm
%
%   ----------------------------------------------------------
%   Acknowledgements: Giovanni Zingaro
%   ----------------------------------------------------------

    % Exact Solution for Height of Brine
    hb_exact = Control.Manhb(Mesh.r, t);

    % Calculate error norms
    Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);
    nq = Quad.nq;
    
    eL2_num = 0;
    eL2_den = 0;
    
    for i = 1:Mesh.ne
        enodes = Mesh.conn(i,:);
        nodesI = Mesh.r(enodes,:);
        
        for p = 1:nq
            N = Quad.Nq{p}';
            dNdxi = Quad.dNdxiq{p};
            wp = Quad.W(p);
            Je = dNdxi'*nodesI;
            

            % quadrature point in physical coordinates
            Ri = nodesI'*N';

            % approximated saturation field at quadrature point
            hb_num_q = N*hb(Mesh.DOF(enodes,:));

            % exact displacement at quadrature point
            hb_exact_q = N*hb_exact(Mesh.DOF(enodes,:));

            % L2 norm
            eL2_num = eL2_num + (hb_num_q - hb_exact_q) * (hb_num_q - hb_exact_q) * wp * Ri * det(Je);
            eL2_den = eL2_den + (hb_exact_q) * (hb_exact_q) * wp * Ri * det(Je);   
           
        end
    end

    m_L2 = sqrt(eL2_num/eL2_den);
    
end


