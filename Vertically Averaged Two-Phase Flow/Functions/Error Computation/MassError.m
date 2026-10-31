function [dhc_scale, MassError] = MassError(Mesh, hc, Qex, Material)


%% Loop through Elements
Qh = 0;
phi = Material.phi;
Sb = Material.sB; % Saturation of Brine [-]

%index = find(hc > 0.0001, 1, 'last');
for e = 1:Mesh.ne
	    
        % Nodal IDs of the Element's Nodes
	    enodes = Mesh.conn(e,:);    
	    
        % Global Coordinates of the Element's Nodes
	    re = Mesh.r(enodes,:); 

        % Elemental Saturation 
        hce = hc(enodes,:);

        % dQ from Saturation over Element
        dQ = (pi*(re(2).^2-re(1).^2))*((hce(2)+hce(1)))*0.5*(phi*(1-Sb));

        % Integrate Saturation of CO2
        Qh = Qh + dQ;
    
end

MassError = abs((Qh-Qex)/Qex);
dhc_scale = 1-(Qh-Qex)/Qex;



