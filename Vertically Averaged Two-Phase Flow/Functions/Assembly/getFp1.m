function Fp1 = getFp1(Mesh, t, BC)
%   GETFP1 -> BOUNDARY FLUX VECTOR
%   Fp1 = getFp1(Mesh, BC, Quad) is a column vector of boundary flux vector 
%   for a axisymmetric vertically averaged multi-phase flow model.
%   
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:   Structure array with the following fields,
%           .nDOF:      Total number of DOFs
%           .r:         Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
% 
%   BC:     Structure array with the following fields,
%           .qc                         Anonymous function of flux boundary
%                                       condition at ends of domain
%           .numActWells                Index of nodes where wells are
%                                       active 


% Acknowledgements: Giovanni Zingaro
  
%% Initialize Boundary Flux Vector
Fp1 = zeros(Mesh.nDOF, 1);

% Loop through Nodes

for i = 1:length(BC.flux_dof)
   
    % Coordinate of Active Well
    Ri = Mesh.r(i,1);

    % DOFn of Active Well
    DOFn = BC.flux_dof(i);
             
    % Injection Flux
    Fp1(DOFn,1) = Fp1(DOFn,1) + (BC.qb(Ri,t)+BC.qc(Ri,t))*Ri;

end
    
end
