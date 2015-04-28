function [G, C, B, I, xyz_to_fd, fd_to_xyz, fd_to_blk, blk_start] = struct_fdd_TSVnew(blks, total_dim, d_value, dz_value, power, t_amb, fluv, TSV)
    
    % discretize the geometry
    total_dim = [total_dim(1)/d_value total_dim(2)/d_value total_dim(3)/dz_value];
    n_blks = size(blks,2);
    TSV_enclosed_area = {};
    %parfor i = 1:n_blks
    for i=1:n_blks
        blks{i}.coordinate = round ([blks{i}.coordinate(1)/d_value, blks{i}.coordinate(2)/d_value, blks{i}.coordinate(3)/dz_value]);
        blks{i}.dim = [blks{i}.dim(1)/d_value, blks{i}.dim(2)/d_value,  blks{i}.dim(3)/dz_value];   
    end
    fd_ind = 0;
    % initiate mapping from xyz index to fd node
    % mapping from fd node to xyz index will be not be initiated since the
    % size is unknown at this time
    
    xyz_to_fd = zeros(round(total_dim(1)), round(total_dim(2)), round(total_dim(3)));
    
    %fluv=0.1620;
    %fluv=1620;
    % shift coordinate, to avoid 0 in index
    parfor i = 1:n_blks
        blks{i}.coordinate = blks{i}.coordinate + [1,1,1];
    end

    % find the corners of the blocks
    x_start = zeros(n_blks,1);
    x_end = zeros(n_blks,1);
    y_start = zeros(n_blks,1);
    y_end = zeros(n_blks,1);
    z_start = zeros(n_blks,1);
    z_end = zeros(n_blks,1);
    blk_start = cell(n_blks,1);
    % blk_end = cell(n_blks,1);
    for i = 1:n_blks
        coordinate = blks{i}.coordinate;
        dim = blks{i}.dim;
        x_start(i) = coordinate(1);
        x_end(i) = coordinate(1)+dim(1)-1;
        y_start(i) = (coordinate(2));
        y_end(i) = (coordinate(2)+dim(2)-1);
        z_start(i) = coordinate(3);
        z_end(i) = coordinate(3)+dim(3)-1;
        blk_start{i} = [x_start(i),y_start(i),z_start(i)];
        % blk_end{i} = [x_end(i),y_end(i),z_end(i)];
    end
    
    % go block by block, build the xyz_to_fd mapping, fd_to_xyz mapping and
    % fd_to_blk mapping
    fd_to_blk = {};
    for i = 1:n_blks

        for l = z_start(i):z_end(i)
            for j = y_start(i):y_end(i)
                for k = x_start(i):x_end(i)
                    % check if the node has already been included in the adjacent
                    % blocks or not
                    is_include = 0;
                    if i > 1
                        for m = 1:i-1
                            if k >= x_start(m) && k <= x_end(m) && j >= y_start(m) && j <= y_end(m)...
                                    && l >= z_start(m) && l <= z_end(m)
                                is_include = 1;
                                % store the block number of the adjacent node
                                fd_ind_tmp = xyz_to_fd(k,j,l);
                                fd_to_blk{fd_ind_tmp} = [fd_to_blk{fd_ind_tmp},i];
                            end
                        end
                    end
                    % if not included by previous blocks, we count this as a
                    % new fd node
                    if is_include == 0
                        fd_ind = fd_ind+1;
                        xyz_to_fd(k,j,l) = fd_ind;
                        fd_to_xyz(fd_ind,:) = [k,j,l];
                        %fprintf('fd ind = %d, x = %d, y = %d, z = %d\n',fd_ind, k,j,l);
                        % store the block number of the node
                        fd_to_blk(fd_ind) = {i};
                    end
                end
            end
        end
    end
    n_fd_node = fd_ind;
    
    % build the matrix based on the mapping
   
    port_ind = 0;
    G = sparse(n_fd_node,n_fd_node);
    C = speye(n_fd_node);   
    for i = 1:n_fd_node

        xyz_ind = fd_to_xyz(i,:);
        
        % Find the adjacent node in xyz indeces, they are also ports
        if xyz_ind(1) == 1 % the most left node does not have left adj
            fd_left = 0;
        else
            % fd_left can be zero or a index value.
            % Zero here means the left one is not in any block, value not set
            % in the xyz_to_fd mapping stage.
            fd_left = xyz_to_fd(xyz_ind(1)-1,xyz_ind(2),xyz_ind(3)); 
            %need to check if fd_left is in an adjacent block
            %If it is in an adjacent block, need to find out how many left
            %grids in adjacent block interface with current grids
            %two more functions are probably needed: (1) fd2xyzMapping (2)
            %find_adjacents
        end
        if xyz_ind(1) == round (total_dim(1)) % most right node
            fd_right = 0;
        else
            fd_right = xyz_to_fd(xyz_ind(1)+1,xyz_ind(2),xyz_ind(3));
        end
        if xyz_ind(2) == 1 % most front node
            fd_front = 0;
        else
            fd_front = xyz_to_fd(xyz_ind(1),xyz_ind(2)-1,xyz_ind(3));
        end
        if xyz_ind(2) == round( total_dim(2)) % most back node
            fd_back = 0;
        else
            fd_back = xyz_to_fd(xyz_ind(1),xyz_ind(2)+1,xyz_ind(3));
        end
        if xyz_ind(3) == 1 % most down node
            fd_down = 0;
        else
            fd_down = xyz_to_fd(xyz_ind(1),xyz_ind(2),xyz_ind(3)-1);
        end
        if xyz_ind(3) == round( total_dim(3)) % most up node
            fd_up = 0;
        else
            fd_up = xyz_to_fd(xyz_ind(1),xyz_ind(2),xyz_ind(3)+1);
        end       
        
        % Stamp the C matrix according to the current node information
        % check the nodes belongs to which block
        fd_cur_blk = fd_to_blk{i};
        rho = blks{fd_cur_blk}.rho;
        cp = blks{fd_cur_blk}.cp;
        C(i,i) = C(i,i)*(rho*cp);
        channel = blks{fd_cur_blk}.channel;
        
                % check if the grid of the node contains TSVs. 
        TSV_locations_inCurrent_Blk = []; 
        if isfield(blks{fd_cur_blk},'TSV_locations') ~= 0
            TSV_locations_inCurrent_Blk = blks{fd_cur_blk}.TSV_locations;
            %nd_left_limit = xyz_ind(1) - 1/2;  nd_right_limit = xyz_ind(1) + 1/2;
            %nd_front_limit = xyz_ind(2) - 1/2;  nd_back_limit = xyz_ind(2) + 1/2;
            nd_left_limit = xyz_ind(1) - 1;  nd_right_limit = xyz_ind(1);
            nd_front_limit = xyz_ind(2);  nd_back_limit = xyz_ind(2) + 1;
 
           % nd_bottom_limit = xyz_ind(3) - 1/2;  nd_top_limit = xyz_ind(3) + 1/2;   
           TSV_cnt = 0;  
           %TSVRadius = TSV.radius * (1-TSV.linear); 
           TSVRadius = TSV.radius;
           %tkddTSV_Cu = TSV.therm_k_Cu/dz_value/dz_value * pi*(TSVRadius)^2/d_value^2 ;
           tkddTSV_Cu = 1/(1/TSV.therm_k_Cu * dz_value/(pi*(TSVRadius)^2)) /(d_value^2*dz_value);
           capTSV_Cu = TSV.rho_Cu*TSV.cp_Cu * pi*(TSVRadius)^2/d_value^2;
%             tkddTSV_Cu = TSV.therm_k_Cu /dz_value /dz_value;
%             capTSV_Cu = TSV.rho_Cu*TSV.cp_Cu ;          
           Vol_Cu = pi * TSVRadius^2  * dz_value;
           %Vol_Ins = pi * TSV.radius^2 * (1 - (1 - TSV.linear)^2) * dz_value;
           Vol_Ins = pi * TSV.radius^2 * ((1 + TSV.linear)^2 - 1) * dz_value;
           
            for loc = 1:length(TSV_locations_inCurrent_Blk)
                if (TSV_locations_inCurrent_Blk{loc}(1) / d_value > nd_left_limit && TSV_locations_inCurrent_Blk{loc}(1) / d_value < nd_right_limit ...
                &&  TSV_locations_inCurrent_Blk{loc}(2) / d_value > nd_front_limit && TSV_locations_inCurrent_Blk{loc}(2) / d_value < nd_back_limit )
               % &&  TSV_locations_inCurrent_Blk{loc}(3) > nd_bottom_limit && TSV_locations_inCurrent_Blk{loc}(2) < nd_top_limit)

%                      %tkddTSV = TSV.therm_k/d_value /d_value; 
%                      TSVRadius = TSV.radius ; 
%                      tkddTSV_Cu = TSV.therm_k_Cu/d_value/d_value * pi*(TSVRadius)^2/d_value^2 ;
%                      capTSV_Cu = TSV.rho_Cu*TSV.cp_Cu * pi*(TSVRadius)^2/d_value^2;
%                      Vol_Cu = pi * TSVRadius^2  * d_value;
%                      Vol_Ins = pi * TSV.radius^2 * ((1 + TSV.linear)^2 - 1) * d_value;
%                     % Vol_bulkSi = d_value^3 - Vol_Cu - Vol_Ins; 
                     TSV_cnt = TSV_cnt + 1;
                     C(i,i) = C(i,i) + (TSV.rho_Cu * TSV.cp_Cu * Vol_Cu + TSV.rho_Ins * TSV.cp_Ins * Vol_Ins )/(d_value^2*dz_value);
                       %  + 0*blks{fd_cur_blk}.rho * blks{fd_cur_blk}.cp * Vol_bulkSi) / d_value^3; 
                     
                    % G(i,i) = G(i,i)+ tkddTSV;
                     % the two parallel resistence in Vertical direction due to
                     % TSVs
                     if fd_up ~=0 
%                          G(i,i) = G(i,i)+ tkddTSV;
%                          G(i,fd_up) = G(i,fd_up) - tkddTSV ; 
%                          C(i,i) = C(i,i) + capTSV;
                         %%% how to adjust G C of silicon in the i-th grid?
                     else
                         %%%step into the empty layer in between until
                         %%%reached the real layer connected by TSV bump.
                         xyz_loc = fd_to_xyz(i,:);
                         xyz_loc(3) = xyz_loc(3) + 1;  %the node above current node, xyz_to_fd should convert this load to 0
                         distance_layers = 1;
                         while(~xyz_to_fd(xyz_loc(1), xyz_loc(2),xyz_loc(3)))
                             xyz_loc(3) = xyz_loc(3) + 1;      
                             distance_layers = distance_layers + 1; 
                         end
                         fd_up_layer = xyz_to_fd(xyz_loc(1), xyz_loc(2),xyz_loc(3));
                         
%                          %%%%% the following lines includes Bump models by using scaling coefficient between cylinder and sphere.  
%                          BumpRadius = sqrt((distance_layers*d_value/2)^2 + TSVRadius^2);
%                          Bumpresist = 1/(TSV.therm_k_Cu*pi*BumpRadius)*log((BumpRadius + sqrt(BumpRadius^2 - TSVRadius^2))/(BumpRadius - sqrt(BumpRadius^2 - TSVRadius^2)));
%                          cylindresist = 1/(TSV.therm_k_Cu)* distance_layers*d_value / (pi*(TSVRadius)^2); 
%                          Resistance_coeff = cylindresist/Bumpresist; 
%                          %BumpV = 2*pi*(BumpRadius^2*(BumpRadius-sqrt(BumpRadius^2-TSVRadius^2)) - (BumpRadius-sqrt(BumpRadius^2-TSVRadius^2))^3/3);
%                          BumpV = 2*pi * (2*BumpRadius^2+TSVRadius^2)/3 * sqrt(BumpRadius^2-TSVRadius^2); 
%                          cylindV = pi*TSVRadius^2*(distance_layers*d_value);
%                          Cap_coeff = BumpV/cylindV; 
%                          G(i,i) = G(i,i)+tkddTSV_Cu / distance_layers * Resistance_coeff;
%                          G(i,fd_up_layer) = G(i,fd_up_layer) - tkddTSV_Cu / distance_layers  * Resistance_coeff;
%                          G(fd_up_layer, i) = G(fd_up_layer,i) - tkddTSV_Cu / distance_layers  * Resistance_coeff;
%                          G(fd_up_layer,fd_up_layer) = G(fd_up_layer, fd_up_layer) + tkddTSV_Cu / distance_layers  * Resistance_coeff;   
%                          C(i,i) = C(i,i) + capTSV_Cu * distance_layers * Cap_coeff ;

%                          % the older method that does not include bump
%                          % shape

                         G(i,i) = G(i,i)+tkddTSV_Cu / distance_layers;
                         G(i,fd_up_layer) = G(i,fd_up_layer) - tkddTSV_Cu / distance_layers ;
                         G(fd_up_layer, i) = G(fd_up_layer,i) - tkddTSV_Cu / distance_layers;
                         G(fd_up_layer,fd_up_layer) = G(fd_up_layer, fd_up_layer) + tkddTSV_Cu / distance_layers;
                                                
                          C(i,i) = C(i,i) + capTSV_Cu * distance_layers ;
                          
%                          tkddInsuLiner =   1/ (-(1/TSV.therm_k_Ins)*log(1 - TSV.linear)/(4*pi*dz_value)) / (d_value^2*dz_value);
%                         G(i,i) = G(i,i)+ tkddInsuLiner;
%                         G(i,fd_up_layer) = G(i,fd_up_layer) - tkddInsuLiner;
%                         G(fd_up_layer, i) = G(fd_up_layer,i) -  tkddInsuLiner;
%                         G(fd_up_layer,fd_up_layer) = G(fd_up_layer, fd_up_layer)  + tkddInsuLiner;


%                          tkddTSV_between_layers = (tkddTSV_Cu / distance_layers)*(tkddInsuLiner)/( tkddTSV_Cu / distance_layers + tkddInsuLiner);
%                          G(i,i) = G(i,i) + tkddTSV_between_layers;
%                          G(i,fd_up_layer) = G(i,fd_up_layer) - tkddTSV_between_layers ;
%                          G(fd_up_layer, i) = G(fd_up_layer,i) - tkddTSV_between_layers;
%                          G(fd_up_layer,fd_up_layer) = G(fd_up_layer, fd_up_layer) + tkddTSV_between_layers;
%                          
                     end
               
                     if fd_down ~=0
%                        G(i,i) = G(i,i)+ tkddTSV;
%                        G(i,fd_down) = G(i,fd_down) - tkddTSV; 
%                        C(i,i) = C(i,i) + capTSV ;
                     end

                end
            end
            %correction of thermal capacitance, the original rho*cp for
            %silicon should be substracted and the portion for silicon in grid i be recalculated
            Vol_bulkSi = d_value^2*dz_value - (Vol_Cu + Vol_Ins)* TSV_cnt; 
            C(i,i) = C(i,i) - (rho*cp) + blks{fd_cur_blk}.rho * blks{fd_cur_blk}.cp * Vol_bulkSi/(d_value^2*dz_value); 
        end

        
        % Stamp the G matrix acording to the resistor with ajacent nodes.
        % the left resistor
        if fd_left ~= 0 && channel == 0
            % Check the block/blocks the left nodes belongs to
            fd_left_blk = fd_to_blk{fd_left};
            r_blks  = [fd_cur_blk,fd_left_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0 
                 tkdd = calculateTSVarrays(r_blks, fd_to_xyz, blks{r_blks(1)}.therm_k, TSV_locations_inCurrent_Blk, blks, i, fd_left, d_value, dz_value, TSV, 'left');
            else 
            AR = blks{r_blks(2)}.dim(3)/blks{r_blks(2)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(2)}.therm_k*Nu; 
            Dh=2*blks{r_blks(2)}.dim(3)*blks{r_blks(2)}.dim(2)/(blks{r_blks(2)}.dim(3)+blks{r_blks(2)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value;
            end
            G(i,i) = G(i,i)+tkdd;
  	        G(i,fd_left) = -tkdd;
        elseif fd_left ~= 0 && channel ~=0
            % Check the block/blocks the left nodes belongs to
            fd_left_blk = fd_to_blk{fd_left};
            r_blks  = [fd_cur_blk,fd_left_blk];
            % calculate the equavalent thermal conductance;
            if blks{r_blks(2)}.channel == 0
            AR = blks{r_blks(1)}.dim(3)/blks{r_blks(1)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(1)}.therm_k*Nu; 
            Dh=2*blks{r_blks(1)}.dim(3)*blks{r_blks(1)}.dim(2)/(blks{r_blks(1)}.dim(3)+blks{r_blks(1)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value;
            else 
            therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
            tkdd = therm_k_eq/d_value/d_value/2;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_left) = -tkdd;
        elseif fd_left ==0 && channel ~=0
            % otherwise, it is inlet of the channel, stamp B and add boundary condition
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
            tkdd = therm_k/d_value/d_value;
            
             G(i,i)=G(i,i)+2*blk.rho*blk.cp*fluv/d_value;
            % G(i,i)=G(i,i)+tkdd+blk.rho*blk.cp*fluv/d_value;
            % boundary current source
            %I(port_ind) = he(5)*d_value/therm_k*tkdd*t_amb;
             I(port_ind) = tkdd*t_amb+2*blk.rho*blk.cp*fluv*t_amb/d_value;
             
        elseif fd_left ==0 && channel ==0
            % otherwise, it is a port, stamp B and add boundary condition
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
%             tkdd = therm_k/d_value/d_value;
%             G(i,i) = G(i,i)+he(5)*d_value/therm_k*tkdd;
            G(i,i) = G(i,i) + he(5) / d_value; 
            % boundary current source
%             I(port_ind) = he(5)*d_value/therm_k*tkdd*t_amb;
            I(port_ind) = he(5)* t_amb / d_value;
        end
        
        % the right resistor
        if fd_right ~= 0 && channel == 0 
            fd_right_blk = fd_to_blk{fd_right};
            r_blks  = [fd_cur_blk,fd_right_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0 
                 tkdd = calculateTSVarrays(r_blks, fd_to_xyz, blks{r_blks(1)}.therm_k, TSV_locations_inCurrent_Blk, blks, i, fd_right, d_value, dz_value, TSV, 'right');
            else 
            AR = blks{r_blks(2)}.dim(3)/blks{r_blks(2)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(2)}.therm_k*Nu; 
            Dh=2*blks{r_blks(2)}.dim(3)*blks{r_blks(2)}.dim(2)/(blks{r_blks(2)}.dim(3)+blks{r_blks(2)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_right) = -tkdd;
        elseif fd_right ~= 0 && channel ~=0
            % Check the block/blocks the left nodes belongs to
            fd_right_blk = fd_to_blk{fd_right};
            r_blks  = [fd_cur_blk,fd_right_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0
            AR = blks{r_blks(1)}.dim(3)/blks{r_blks(1)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(1)}.therm_k*Nu; 
            Dh=2*blks{r_blks(1)}.dim(3)*blks{r_blks(1)}.dim(2)/(blks{r_blks(1)}.dim(3)+blks{r_blks(1)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value;               
            else 
            therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
            tkdd = therm_k_eq/d_value/d_value/2;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_right) = -tkdd;
        elseif fd_right ==0 && channel ~=0
            % otherwise, it is the outlet of the channel 
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
            tkdd = therm_k/d_value/d_value;
            % if using Neuman boundary condition at outlet
            %G(i,i)=G(i,i)-tkdd;
            G(i,i)=G(i,i)-tkdd + blk.rho*blk.cp*fluv/d_value;
            G(i,fd_left)=-tkdd-blk.rho*blk.cp*fluv/d_value;
            % boundary current source
            % if using Neuman boundary condition at outlet
            I(port_ind)=tkdd*t_amb*0;
         elseif fd_right == 0 && channel == 0
            % otherwise, it is a port
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
%             tkdd = therm_k/d_value/d_value;
%             G(i,i) = G(i,i)+he(6)*d_value/therm_k*tkdd;
            G(i,i) = G(i,i) + he(6) / d_value; 
            % boundary current source
%             I(port_ind) = he(6)*d_value/therm_k*tkdd*t_amb;
            I(port_ind) = he(6)* t_amb / d_value;
        end
        
        % the front resistor
        if fd_front ~= 0 && channel ==0
            fd_front_blk = fd_to_blk{fd_front};
            r_blks  = [fd_cur_blk,fd_front_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0   
                 tkdd = calculateTSVarrays(r_blks, fd_to_xyz, blks{r_blks(1)}.therm_k, TSV_locations_inCurrent_Blk, blks, i, fd_front, d_value, dz_value,TSV, 'front');
            else 
            AR = blks{r_blks(2)}.dim(3)/blks{r_blks(2)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(2)}.therm_k*Nu; 
            Dh=2*blks{r_blks(2)}.dim(3)*blks{r_blks(2)}.dim(2)/(blks{r_blks(2)}.dim(3)+blks{r_blks(2)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_front) = -tkdd;
        elseif fd_front ~= 0 && channel ~= 0
            fd_front_blk = fd_to_blk{fd_front};
            r_blks  = [fd_cur_blk,fd_front_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0
            AR = blks{r_blks(1)}.dim(3)/blks{r_blks(1)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(1)}.therm_k*Nu; 
            Dh=2*blks{r_blks(1)}.dim(3)*blks{r_blks(1)}.dim(2)/(blks{r_blks(1)}.dim(3)+blks{r_blks(1)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value; 
            else 
            therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
            tkdd = therm_k_eq/d_value/d_value/2;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_front) = -tkdd;
        elseif fd_front == 0 % channel cannot be at the boundary
            % otherwise, it is a port
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
%             tkdd = therm_k/d_value/d_value;
%             G(i,i) = G(i,i)+he(3)*d_value/therm_k*tkdd;
            G(i,i) = G(i,i) + he(3) / d_value; 
            % boundary current source
%             I(port_ind) = he(3)*d_value/therm_k*tkdd*t_amb;
            I(port_ind) = he(3)* t_amb / d_value;
        end
        
        % the back resistor
        if fd_back ~= 0 && channel == 0
            fd_back_blk = fd_to_blk{fd_back};
            % same_elem contains the blocks # which the resistor belongs to
            r_blks  = [fd_cur_blk,fd_back_blk];
            % calculate the equavalent thermal conductance
            AR = blks{r_blks(2)}.dim(3)/blks{r_blks(2)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            if blks{r_blks(2)}.channel == 0   
                 tkdd = calculateTSVarrays(r_blks, fd_to_xyz, blks{r_blks(1)}.therm_k, TSV_locations_inCurrent_Blk, blks, i, fd_back, d_value, dz_value, TSV, 'back');
            else 
            AR = blks{r_blks(2)}.dim(3)/blks{r_blks(2)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(2)}.therm_k*Nu; 
            Dh=2*blks{r_blks(2)}.dim(3)*blks{r_blks(2)}.dim(2)/(blks{r_blks(2)}.dim(3)+blks{r_blks(2)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_back) = -tkdd;
       elseif fd_back ~= 0 && channel ~= 0
            fd_back_blk = fd_to_blk{fd_back};
            r_blks  = [fd_cur_blk,fd_back_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0
            AR = blks{r_blks(1)}.dim(3)/blks{r_blks(1)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(1)}.therm_k*Nu; 
            Dh=2*blks{r_blks(1)}.dim(3)*blks{r_blks(1)}.dim(2)/(blks{r_blks(1)}.dim(3)+blks{r_blks(1)}.dim(2))*d_value;
            tkdd = therm_k_eq/Dh/d_value;       
            else 
            therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
            tkdd = therm_k_eq/d_value/d_value/2;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_back) = -tkdd;
        elseif fd_back == 0
            % otherwise, it is a port
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
%             tkdd = therm_k/d_value/d_value;
%             G(i,i) = G(i,i)+he(4)*d_value/therm_k*tkdd;
            G(i,i) = G(i,i) + he(4) / d_value; 
            % boundary current source
%             I(port_ind) = he(4)*d_value/therm_k*tkdd*t_amb;
            I(port_ind) = he(4)* t_amb / d_value;
        end
        
        % the up resistor
        if fd_up ~= 0 && channel == 0
            fd_up_blk = fd_to_blk{fd_up};
            % same_elem contains the blocks # which the resistor belongs to
            r_blks  = [fd_cur_blk,fd_up_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0 
                 tkdd = calculateTSV_Vpath(r_blks, fd_to_xyz, blks{r_blks(1)}.therm_k, TSV_locations_inCurrent_Blk, blks, i, fd_up, d_value, dz_value, TSV, 'up');
            else 
            AR = blks{r_blks(2)}.dim(3)/blks{r_blks(2)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(2)}.therm_k*Nu; 
            Dh=2*blks{r_blks(2)}.dim(3)*blks{r_blks(2)}.dim(2)/(blks{r_blks(2)}.dim(3)+blks{r_blks(2)}.dim(2))*dz_value;
            tkdd = therm_k_eq/Dh/dz_value;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_up) = -tkdd;
       elseif fd_up ~= 0 && channel ~= 0
            fd_up_blk = fd_to_blk{fd_up};
            r_blks  = [fd_cur_blk,fd_up_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0
            AR = blks{r_blks(1)}.dim(3)/blks{r_blks(1)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(1)}.therm_k*Nu; 
            Dh=2*blks{r_blks(1)}.dim(3)*blks{r_blks(1)}.dim(2)/(blks{r_blks(1)}.dim(3)+blks{r_blks(1)}.dim(2))*dz_value;
            tkdd = therm_k_eq/Dh/dz_value;       
            else 
             therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
             tkdd = therm_k_eq/dz_value/dz_value/2;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_up) = -tkdd;
        elseif fd_up == 0
            % otherwise, it is a port
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
%             tkdd = therm_k/d_value/d_value; 
%             G(i,i) = G(i,i)+he(2)*d_value/therm_k*tkdd;
            G(i,i) = G(i,i) + he(2) / dz_value; 
            % boundary current source
%             I(port_ind) = he(2)*d_value/therm_k*tkdd*t_amb;
            I(port_ind) = he(2)* t_amb / dz_value;
        end
        
        % the down resistor
        if fd_down ~= 0 && channel == 0
            fd_down_blk = fd_to_blk{fd_down};
            % same_elem contains the blocks # which the resistor belongs to
            r_blks  = [fd_cur_blk,fd_down_blk];
            % calculate the equavalent thermal conductance;
            if blks{r_blks(2)}.channel == 0
                 tkdd = calculateTSV_Vpath(r_blks, fd_to_xyz, blks{r_blks(1)}.therm_k, TSV_locations_inCurrent_Blk, blks, i, fd_down, d_value, dz_value, TSV, 'down');
            else 
            AR = blks{r_blks(2)}.dim(3)/blks{r_blks(2)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(2)}.therm_k*Nu; 
            Dh=2*blks{r_blks(2)}.dim(3)*blks{r_blks(2)}.dim(2)/(blks{r_blks(2)}.dim(3)+blks{r_blks(2)}.dim(2))*dz_value;
            tkdd = therm_k_eq/Dh/dz_value;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_down) = -tkdd;
       elseif fd_down ~= 0 && channel ~= 0
            fd_down_blk = fd_to_blk{fd_down};
            r_blks  = [fd_cur_blk,fd_down_blk];
            % calculate the equavalent thermal conductance
            if blks{r_blks(2)}.channel == 0
            AR = blks{r_blks(1)}.dim(3)/blks{r_blks(1)}.dim(2);
            Nu = 8.235*(1-2.0421*AR+3.0853*AR^2-2.4765*AR^3+1.0578*AR^4-0.1861*AR^5);
            therm_k_eq =  blks{r_blks(1)}.therm_k*Nu; 
            Dh=2*blks{r_blks(1)}.dim(3)*blks{r_blks(1)}.dim(2)/(blks{r_blks(1)}.dim(3)+blks{r_blks(1)}.dim(2))*dz_value;
            tkdd = therm_k_eq/Dh/dz_value;       
            else 
            therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
            tkdd = therm_k_eq/dz_value/dz_value/2;
            end
            G(i,i) = G(i,i)+tkdd;
            G(i,fd_down) = -tkdd;
        elseif fd_down == 0
            % otherwise, it is a port
            port_ind = port_ind+1;
            b_row(port_ind) = i;
            % boundary resistor to GND
            blk = blks{fd_to_blk{i}};
            therm_k = blk.therm_k;
            he = blk.he;
%             tkdd = therm_k/d_value/d_value;
%             G(i,i) = G(i,i)+he(1)*d_value/therm_k*tkdd;
            G(i,i) = G(i,i) + he(1) / dz_value; 
            % boundary current source
%            I(port_ind) = he(1)*d_value/therm_k*tkdd*t_amb;
            I(port_ind) = he(1)* t_amb / dz_value;
        end
         blk=blks{fd_to_blk{i}};
         %v=1.620;
         C_conv=fluv*blk.rho*blk.cp/d_value/2;
        if channel ~= 0 && fd_right ~=0 && fd_left ~=0
         % + - defines the flow of the fluidics
         % the fluidics flows from inlet to outlet and its not reversable given the boundary condition
         G(i,fd_right)=G(i,fd_right)+C_conv;
	     G(i,fd_left)=G(i,fd_left)-C_conv;
        end

         
    end

    power_idx = power.idx;
    power_val = power.value;
    n_power = length(power_val);
    for i = 1:n_power
        power_fd = xyz_to_fd(power_idx(i,1),power_idx(i,2),power_idx(i,3));
        if power_fd == 0
            error('power source is out of the block');
        end
        port_ind = port_ind+1;
        b_row(port_ind) = power_fd;
        I(port_ind) = power_val(i);
    end
    
    n_port = port_ind;
    B = sparse(n_fd_node,n_port);
    for i = 1:n_port
        B(b_row(i),i) = 1;
    end 

% function same_elem = find_same_elem(array_1,array_2);
%     leng_1 = length(array_1);
%     leng_2 = length(array_2);
%     same_elem = [];
%     for i = 1:leng_1
%         for j = 1:leng_2
%             if array_1(i) == array_2(j)
%                 same_elem = [same_elem array_1(i)];
%                 break;
%             end
%         end
%     end
