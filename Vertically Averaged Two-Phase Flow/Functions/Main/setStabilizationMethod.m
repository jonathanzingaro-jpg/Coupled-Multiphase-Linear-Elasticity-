function [storagestabfile_name, diffusionstabfile_name] = setStabilizationMethod(Control)
%   setStabilizationMethod is used to initialize pointer function based on
%   the stabilization method. 
%   SU -> Control.Stab = 1, LSFEM = Control.Stab = 2
%   ----------------------------------------------------------------------
%   Created by: Jonathan Zingaro
%   ----------------------------------------------------------------------

    switch Control.Stab
        case 1
            storagestabfile_name = 'getCsu';
            diffusionstabfile_name = 'getKsu';
        case 2
            storagestabfile_name = 'getCgls';
            diffusionstabfile_name = 'getKgls';
        case 3
            storagestabfile_name = 'getCsupg';
            diffusionstabfile_name = 'getKsupg';
    end

end

