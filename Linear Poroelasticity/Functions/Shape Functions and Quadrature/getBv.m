function Bv = getBv(B,N,r,ElementType)
%GETBV Shape function derivative matrix in Voigt form
%   Bv = GETBV(B, dim) transforms the matrix of shape function 
%   derivatives, B, to Voigt notation, returning a matrix with n*dim rows 
%   in which n is the number of element nodes and dim is the number of 
%   spatial dimensions). It should be mentioned this constructs the B
%   matrix for axisymmetric elements. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   B:      vector of shape function derivatives (size n x dim where n 
%           is the number of element nodes)
%   dim:    number of spatial dimensions
% 
%   --------------------------------------------------------------------
%   Example
%   --------------------------------------------------------------------
%   B = |N1,1 N2,1 ... Nn,1|^T
%       |N1,2 N2,2 ... Nn,2|

%   Bv = |N1,1 0     N2,1 0     ...  Nn,1 0   |^T
%        |0    N1,2  0    N2,2  ...  0    Nn,2|
%        |N1,2 N1,1  N2,2 N2,1  ...  Nn,2 Nn,1|

%   Acknowledgements: Matin Parchei Esfahani, Jack Chessa, Jonathan Zingaro

    % Construct B Matrix in Voigt Notation for Axisymmetric Element for Q4
    % Elements
    
    switch ElementType
        case 'Q4'
            Bv = [B(1,1),0,B(1,2),0,B(1,3),0,B(1,4),0; ...
                  0,B(2,1),0,B(2,2),0,B(2,3),0,B(2,4); ...
                  N(1)/r,0,N(2)/r,0,N(3)/r,0,N(4)/r,0; ...
                  B(2,1),B(1,1),B(2,2),B(1,2),B(2,3),B(1,3),B(2,4),B(1,4)];
        case 'T3'
             Bv = [B(1,1),0,B(1,2),0,B(1,3),0; ...
                   0,B(2,1),0,B(2,2),0,B(2,3); ...
                   N(1)/r,0,N(2)/r,0,N(3)/r,0; ...
                   B(2,1),B(1,1),B(2,2),B(1,2),B(2,3),B(1,3)];
        case 'T6'
            Bv = [B(1,1), 0, B(1,2), 0, B(1,3), 0, B(1,4), 0, B(1,4), 0, B(1,4),0; ...
                  0, B(2,1), 0, B(2,2), 0, B(2,3), 0, B(2,4), 0, B(2,5), 0, B(2,6); ...
                  N(1)/r, 0, N(2)/r, 0, N(3)/r, 0, N(4)/r, 0, N(5)/r, 0, N(6)/r, 0; ...
                  B(2,1), B(1,1), B(2,2), B(1,2), B(2,3), B(1,3), B(2,4), B(1,4), B(2,5), B(1,5), B(2,6), B(1,6)];
        case 'Q9'
             Bv = [B(1,1), 0, B(1,2), 0, B(1,3), 0, B(1,4), 0, B(1,5), 0, B(1,6), 0, B(1,7), 0, B(1,8), 0, B(1,9), 0; ...
                   0, B(2,1), 0, B(2,2), 0, B(2,3), 0, B(2,4), 0, B(2,5), 0, B(2,6), 0, B(2,7), 0, B(2,8), 0, B(2,9); ...
                   N(1)/r, 0, N(2)/r, 0, N(3)/r, 0, N(4)/r, 0, N(5)/r, 0, N(6)/r, 0, N(7)/r, 0, N(8)/r, 0, N(9)/r, 0; ...
                   B(2,1), B(1,1), B(2,2), B(1,2), B(2,3), B(1,3), B(2,4), B(1,4), B(2,5), B(1,5), B(2,6), B(1,6), B(2,7), B(1,7), B(2,8), B(1,8), B(2,9), B(1,9)]; 
        case 'Q8'
             Bv = [B(1,1), 0, B(1,2), 0, B(1,3), 0, B(1,4), 0, B(1,5), 0, B(1,6), 0, B(1,7), 0, B(1,8), 0; ...
                   0, B(2,1), 0, B(2,2), 0, B(2,3), 0, B(2,4), 0, B(2,5), 0, B(2,6), 0, B(2,7), 0, B(2,8); ...
                   N(1)/r, 0, N(2)/r, 0, N(3)/r, 0, N(4)/r, 0, N(5)/r, 0, N(6)/r, 0, N(7)/r, 0, N(8)/r, 0; ...
                   B(2,1), B(1,1), B(2,2), B(1,2), B(2,3), B(1,3), B(2,4), B(1,4), B(2,5), B(1,5), B(2,6), B(1,6), B(2,7), B(1,7), B(2,8), B(1,8)]; 
    end
end
