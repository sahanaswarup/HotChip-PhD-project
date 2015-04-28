function [glob_dat, blk_dat] = read_file_TSV(file_name)

    fid = fopen(file_name,'r');
    if fid == -1
        error('Cannot open file %s! Please check the file path and permission.', file_name);
    end
    
    blk_idx = 0;
    blk_flag = 0; % flag for block region
    is_in_comment = 0; % flag for comment region
    while 1
        tline = fgetl(fid);
        if tline == -1
            break;
        end
        [token, remain] = strtok(tline);
        % skip the empty line
        if length(tline) == 0 
            continue;
        end
        % skip the # comment
        if strcmp(token(1),'#')
            continue;
        end
        % skip the /* */ comment
        if strcmp(token(1),'/') && strcmp(token(2),'*')
            if is_in_comment == 0
                is_in_comment = 1; % begin comment
                continue;
            else
                error('/* used after /*, consider using */ to end comment');
            end
        end
        if strcmp(token(1),'*') && strcmp(token(2),'/')
            if is_in_comment == 1
                is_in_comment = 0; % end comment
                continue;
            else
                error('*/ used before /*, consider using /* to begin comment');
            end
        end
        if is_in_comment == 1 % between /* and */, skip
            continue;
        end
        
        if blk_flag == 0 % parsing global info, until 'block_begin'
            switch token
              case {'scale'}
                [token, remain] = strtok(remain);
                glob_dat.scale = str2double(token);
              case {'d_value'}
                [token, remain] = strtok(remain);
                glob_dat.d_value = str2double(token);
              case {'dz_value'}
                [token, remain] = strtok(remain);
                glob_dat.dz_value = str2double(token);
              case {'t_amb'}
                [token, remain] = strtok(remain);
                glob_dat.t_amb = str2double(token);
              case {'TSV_therm_k_Cu'}
                [token, remain] = strtok(remain);
                glob_dat.TSV_therm_k_Cu = str2double(token);
              case {'TSV_rho_Cu'}
                [token, remain] = strtok(remain);
                glob_dat.TSV_rho_Cu = str2double(token); 
              case {'TSV_cp_Cu'} 
                [token, remain] = strtok(remain);
                glob_dat.TSV_cp_Cu = str2double(token); 
              case {'TSV_therm_k_Ins'}
                [token, remain] = strtok(remain);
                glob_dat.TSV_therm_k_Ins = str2double(token);
              case {'TSV_rho_Ins'}
                [token, remain] = strtok(remain);
                glob_dat.TSV_rho_Ins = str2double(token); 
              case {'TSV_cp_Ins'} 
                [token, remain] = strtok(remain);
                glob_dat.TSV_cp_Ins = str2double(token);
              case {'TSV_radius'} 
                [token, remain] = strtok(remain);
                glob_dat.TSV_radius = str2double(token); 
              case {'TSV_linear'} 
                [token, remain] = strtok(remain);
                glob_dat.TSV_linear = str2double(token); 
              case {'block_begin'} % blk detected, change to blk parsing
                blk_idx = blk_idx+1;
                blk_flag = 1;
            end
        else 
            % parsing blk info, until 'block_end'
            switch token
              case {'coordinate'}
                tmp = textscan(remain,'%f');
                blk_dat{blk_idx}.coordinate = tmp{1}';
              case {'dim'}
                tmp = textscan(remain,'%f');
                blk_dat{blk_idx}.dim = tmp{1}';
              case {'therm_k'}
                blk_dat{blk_idx}.therm_k = str2double(remain);
              case {'he'}
                tmp = textscan(remain,'%f');
                blk_dat{blk_idx}.he = tmp{1}';
              case {'rho'}
                blk_dat{blk_idx}.rho = str2double(remain);
              case {'cp'}
                blk_dat{blk_idx}.cp = str2double(remain);
              case {'channel'}
                blk_dat{blk_idx}.channel =  str2double(remain);
              case {'TSV_locations'}
                tmp = textscan(remain,'%f');
                if isfield(blk_dat{blk_idx},'TSV_locations') == 0
                   blk_dat{blk_idx}.TSV_locations = {tmp{1}'};
                else
                   blk_dat{blk_idx}.TSV_locations = [blk_dat{blk_idx}.TSV_locations, tmp{1}'];
                end
              case {'block_end'}
                % finished parsing the current blk, change to global parsing
                blk_flag = 0;
            end
        end
    end
    fclose(fid);
        
    % check if any error, and get total_dim information
    [total_dim, blk_dat] = struct_error_check(glob_dat, blk_dat);
    glob_dat.total_dim = total_dim;
    
function [total_dim, blks] = struct_error_check(glob_dat, blks)
    
    d_value = glob_dat.d_value;
    n_blks = length(blks);
    % try discretization here
%     for i = 1:n_blks
%         blk_coordinate_tmp = blks{i}.blk_coordinate./d_value;
%         dim_tmp = blks{i}.dim./d_value;
%         %if nnz(mod(blk_coordinate_tmp,1))~=0 || nnz(mod(dim_tmp,1))~=0
%         if nnz(mod(dim_tmp,1))~=0
%             error('The blocks indeces cannot be dircretized into integers, please change d_value');
%         end        
%     end
    
    % check the overlapping of the blocks       
    x_start = zeros(n_blks,1);
    x_end = zeros(n_blks,1);
    y_start = zeros(n_blks,1);
    y_end = zeros(n_blks,1);
    z_start = zeros(n_blks,1);
    z_end = zeros(n_blks,1);
    % now, check from block to block, corner to corner
    for i = 1:n_blks
        blk_coordinate = blks{i}.coordinate;
        blk_coordinate_x(i) = blk_coordinate(1);
        blk_coordinate_y(i) = blk_coordinate(2);
        blk_coordinate_z(i) = blk_coordinate(3);
        dim = blks{i}.dim;
        x_start(i) = blk_coordinate(1);
        x_end(i) = blk_coordinate(1) + dim(1);
        y_start(i) = blk_coordinate(2);
        y_end(i) = blk_coordinate(2)+dim(2);
        z_start(i) = blk_coordinate(3);
        z_end(i) = blk_coordinate(3)+dim(3);
        % compare i-th block with all the previous blocks
        if i > 1
            for j = 1:i-1 
                if blk_coordinate_x(i) >= blk_coordinate_x(j) && blk_coordinate_y(i) >= blk_coordinate_y(j) && blk_coordinate_z(i) >= blk_coordinate_z(j)
                    if x_start(i) < x_end(j) && y_start(i) < y_end(j) && z_start(i) < z_end(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end
                if blk_coordinate_x(i) >= blk_coordinate_x(j) && blk_coordinate_y(i) >= blk_coordinate_y(j) && blk_coordinate_z(i) <= blk_coordinate_z(j)
                    if x_start(i) < x_end(j) && y_start(i) < y_end(j) && z_end(i) > z_start(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end
                if blk_coordinate_x(i) >= blk_coordinate_x(j) && blk_coordinate_y(i) <= blk_coordinate_y(j) && blk_coordinate_z(i) >= blk_coordinate_z(j)
                    if x_start(i) < x_end(j) && y_end(i) > y_start(j) && z_start(i) < z_end(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end
                if blk_coordinate_x(i) <= blk_coordinate_x(j) && blk_coordinate_y(i) >= blk_coordinate_y(j) && blk_coordinate_z(i) >= blk_coordinate_z(j)
                    if x_end(i) > x_start(j) && y_start(i) < y_end(j) && z_start(i) < z_end(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end        
                if blk_coordinate_x(i) >= blk_coordinate_x(j) && blk_coordinate_y(i) <= blk_coordinate_y(j) && blk_coordinate_z(i) <= blk_coordinate_z(j)
                    if x_start(i) < x_end(j) && y_end(i) > y_start(j) && z_end(i) > z_start(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end
                if blk_coordinate_x(i) <= blk_coordinate_x(j) && blk_coordinate_y(i) >= blk_coordinate_y(j) && blk_coordinate_z(i) <= blk_coordinate_z(j)
                    if x_end(i) > x_start(j) && y_start(i) < y_end(j) && z_end(i) > z_start(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end
                if blk_coordinate_x(i) <= blk_coordinate_x(j) && blk_coordinate_y(i) <= blk_coordinate_y(j) && blk_coordinate_z(i) >= blk_coordinate_z(j)
                    if x_end(i) > x_start(j) && y_end(i) > y_start(j) && z_start(i) < z_end(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end
                if blk_coordinate_x(i) <= blk_coordinate_x(j) && blk_coordinate_y(i) <= blk_coordinate_y(j) && blk_coordinate_z(i) <= blk_coordinate_z(j)
                    if x_end(i) > x_start(j) && y_end(i) > y_start(j) && z_end(i) > z_start(j)
                        error('block %d inside block %d, please check your str file', i,j);
                    end
                end
            end
        end
    end
    
    % check if any block out of the first quadrant and calculate
    % total_dim for output
    [x_start_min,i] = min(x_start);
    if x_start_min < 0
        error('block %d out of the first quadrant', i);
    end
    x_end_max = max(x_end);
    [y_start_min,i] = min(y_start);
    if y_start_min < 0
        error('block %d out of the first quadrant', i);
    end
    y_end_max = max(y_end);
    [z_start_min,i] = min(z_start);
    if z_start_min < 0
        error('block %d out of the first quadrant', i);
    end
    z_end_max = max(z_end);
    
    total_dim = [x_end_max-x_start_min, y_end_max-y_start_min, z_end_max-z_start_min];
    
    offset = [x_start_min, y_start_min, z_start_min];

    for i = 1:n_blks
        blks{i}.coordinate = blks{i}.coordinate-offset;
    end
