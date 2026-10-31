function writeData2txt_Disp(d, Mesh, data_dir, BC)
%writeData2txt Exports results to data as a textfile located in the
%   data_dir
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:       Structure array with the following fields, 
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .conn:  Array of element connectivity (size ne x nne)
%               .xdofs: DOFs in the x-direction
%               .ydofs: DOFs in the y-direction
%               .zdofs: DOFs in the z-direction
%               .nn:    Total number of nodes 
%               .nsd:   Number of spatial dimensions
%               .nDOF:  Total number of DOFs
%               .ne:    Total number of elements in the mesh
%   d:          Column vector of displacements (size ndof x 1 where ndof 
%               is the number of degrees of freedom)
%   Pressure:   Structure containing pressure parameters
%   Material:   Structure containing material parameters
%   BC:         Structure containing node sets
%
%   Acknowledgement: Jonathan Zingaro

    % Header information
    header = "Displacement at Measurement Points [mm]";

    % Convert displacement values to strings with specified format
    dm_values = d(BC.measurement) .* 1000;
    dm_str = compose("%.16E", dm_values);
    
    % Combine header and formatted displacement strings
    dm_data = [header; dm_str];
    
    % Writing the data to a text file
    writematrix(dm_data, fullfile(data_dir, 'dm.txt'));

    

end
