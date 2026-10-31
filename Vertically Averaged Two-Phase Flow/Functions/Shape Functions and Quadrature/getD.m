function D = getD(type)
%GETD calculates the second derivative of the shape functions used in SUPG
%   and LSFEM for the saturation equation.
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   D:      vector of shape function second derivatives (size n x dim 
%           where n is the number of element nodes)
%   Mesh.type: type of elements used in the spatial approximation
%   --------------------------------------------------------------------
%   Acknowledgements: Jonathan Zingaro

switch type
    case 'L2'
        D = [0;0];
    case 'L3'
        D = [1; -2; 1];

end
end
