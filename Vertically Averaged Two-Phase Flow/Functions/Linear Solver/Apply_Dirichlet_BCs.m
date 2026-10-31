function [F,K]=Apply_Dirichlet_BCs(fixunode,fixuvalue,F,K, progress_on)

Ko=K;
% Apply BCs
if progress_on == 1
    disp([num2str(toc),': Applying Dirichlet BCs...']);
else
end
nfix = length(fixunode); % number of nodes where u BCs will be applied.
nfixeddof = 1*nfix; % number of degrees of freedom to be prescribed - 1DOF per node.
% fixed_u_dofs = zeros(nfixeddof,1); 
fixed_u_dofs=fixunode;

if isempty(fixed_u_dofs) == 0
    for ii=1:length(fixed_u_dofs)
        F=F-Ko(:,fixed_u_dofs(ii))*fixuvalue(ii);
        F(fixed_u_dofs(ii))=fixuvalue(ii);
        
        % zero out the rows and columns of the K matrix
        K(fixed_u_dofs(ii),:)=0;
        K(:,fixed_u_dofs(ii))=0;
        K(fixed_u_dofs(ii),fixed_u_dofs(ii))=speye(length(fixed_u_dofs(ii)));
    end
end

end

