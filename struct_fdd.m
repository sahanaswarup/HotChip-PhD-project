function [heat_map, min_temp, max_temp] = struct_fdd(filename, power_loc_file)
%% read in the configuration
%[glob_dat, blks] = read_file_TSV('./chip_struct_TSV3');
%[glob_dat, blks] = read_file_TSV('./example_TSV_newmodelscaled3.str');
%[glob_dat, blks] = read_file_TSV('./example_TSV_newmodel_paper.str');
[glob_dat, blks] = read_file_TSV(filename); % equivalent to what's in the GUI
scale = glob_dat.scale;
total_dim = glob_dat.total_dim./scale;
d_value = glob_dat.d_value./scale;
dz_value = glob_dat.dz_value./scale;
TSV_radius = glob_dat.TSV_radius./scale;
TSV_linear = glob_dat.TSV_linear; %it is a ratio
t_amb = 273.15+glob_dat.t_amb;
n_blk = length(blks);
for i = 1:n_blk
    blks{i}.coordinate = blks{i}.coordinate./scale;
    blks{i}.dim = blks{i}.dim./scale;
    if isfield(blks{i},'TSV_locations') ~= 0
        for loc = 1:length(blks{i}.TSV_locations)
          blks{i}.TSV_locations{loc} = blks{i}.TSV_locations{loc} ./ scale;
        end
    end
end
%TSV = struct('therm_k', glob_dat.TSV_therm_k, 'rho',glob_dat.TSV_rho, 'cp', glob_dat.TSV_cp, 'cross_area', TSV_cross_area); 
TSV = struct('therm_k_Cu', glob_dat.TSV_therm_k_Cu, 'rho_Cu',glob_dat.TSV_rho_Cu, 'cp_Cu', glob_dat.TSV_cp_Cu, ...
             'therm_k_Ins', glob_dat.TSV_therm_k_Ins, 'rho_Ins', glob_dat.TSV_rho_Ins, 'cp_Ins', glob_dat.TSV_cp_Ins, ...
             'radius', TSV_radius, 'linear', TSV_linear); 
%% put power here

total_blks1 = blks{1}.dim./d_value;
% If first layer coordinate is not (0,0,0), find offset. Assumes that if
% layers are not at offset, then heat sink is. If heat sink and layers are
% all starting at offset, then this logic is screwed. Ideally, we need to
% find the minimum x,y,z coordinate across all layers and then take offset
% of coordinate of first block from that min and add it to the power pos.
% Power pos specifies the block id within the layer only. This way it
% doesnt change depending on whether heat sink is larger or same as the
% layers
offset      = blks{1}.coordinate./d_value;
%x_p1 = ceil(total_blks1(1)/4);
%y_p1 = ceil(total_blks1(2)/4);
%z_p1 = ceil(total_blks1(3)/2);
%x_p2 = ceil(total_blks1(1)*1/4);
%y_p2 = ceil(total_blks1(2)*1/4);
%z_p2 = ceil(total_blks1(3)/2)+2;
%x_p3 = ceil(total_blks1(1)/4)+1;
%y_p3 = ceil(total_blks1(2)/4)+1;
%z_p3 = ceil(total_blks1(3)/2);
%x_p4 = ceil(total_blks1(1)*3/4)+1;
%y_p4 = ceil(total_blks1(2)*3/4)+1;
%z_p4 = ceil(total_blks1(3)/2);

% power_idx = [x_p1,y_p1,z_p1;
%     %    x_p2,y_p2,z_p2;
%         x_p3,y_p3,z_p3;
%         x_p4,y_p4,z_p4];
%z_p = ceil(total_blks1(3)/2);
%z_p = 1;
%i=2;
%for i = 1:total_blks1(2)
%        %pp = i+(j-1)*total_blks1(1);%+(z_p-1)*total_blks1(1)*total_blks1(2);
%        power_idx(j,:) = [i, j, z_p];
%end;

power_fp = csvread(power_loc_file);
z_p = 1;
for i = 1:size(power_fp, 1)
    x_pos = power_fp(i,1);
    y_pos = power_fp(i,2);
    % The index of the power index has to be w.r.t heat sink origin. If
    % heat sink is larger than the layers, then fd_to_xyz is allocated with
    % the larger size. Therefore the x_pos and y_pos need to be specified
    % properly. We need to think how this can be simplified
    
    if x_pos <= 0 || x_pos > total_blks1(1)
        continue
    end
    if y_pos <= 0 || y_pos > total_blks1(2)
        continue  
    end
    power_idx(i, :) = [x_pos + offset(1), y_pos + offset(1), z_p];
end
% j=2;
% for i = 1:total_blks1(1)
%         pp = i+(j-1)*total_blks1(1);%+(z_p-1)*total_blks1(1)*total_blks1(2);
%         power_idx(pp,:) = [i, j, z_p];
% end;


% i=3;
% for j = 1:total_blks1(2)
%         pp = i+(j-1)*total_blks1(1);%+(z_p-1)*total_blks1(1)*total_blks1(2);
%         power_idx(pp,:) = [i, j, z_p];
% end;
% pp_tmp = pp;
% %%%used for 2-layer silcon
% for j = 1:total_blks1(2)
%         pp = pp_tmp+i+(j-1)*total_blks1(1);
%         power_idx(pp,:) = [i, j, 14+z_p];
% end;
%remove zeros
%p1=power_idx((power_idx(:,1) ~=0),1);
%p2=power_idx((power_idx(:,2) ~=0),2);
%p3=power_idx((power_idx(:,3) ~=0),3);
%power_idx = [p1, p2, p3];

%power_val = 3.2e11/8/4;
% power_val = [power_val,power_val,power_val];
% power = struct('idx',power_idx,'value',power_val);
% n_power = length(power_val);

power_den2D = 5e6; %5e6;
power_val = power_den2D/dz_value;  %power_density is W/m^2, need to divide d_value to get to W/m^3
p_density=power_val;
power_val = repmat(power_val,1,size(power_idx,1));
power = struct('idx',power_idx,'value',power_val);
n_power = length(power_val);

%% generate the fdd model
% IMPORTANT - equivalent to gui - builds gcb matrix
disp('Building thermal model...');
v_flow = 1620;
tic;
[G, C, B, I, xyz_to_fd, fd_to_xyz, fd_to_blk, blk_start] = struct_fdd_TSVnew(blks, total_dim, d_value, dz_value, power, t_amb, v_flow, TSV);

%[G, C, B, I, xyz_to_fd, fd_to_xyz, fd_to_blk, blk_start] = struct_fdd_paper3bak(blks, total_dim, d_value, power, t_amb, v_flow);

struct_fdd_time = toc
% % display the structure of the state space matrices
% figure;
% spy(G);
% figure;
% spy(C);
% figure;
% spy(B);
x_position = fd_to_xyz(:,1);
y_position = fd_to_xyz(:,2);
z_position = fd_to_xyz(:,3);
%plot3(x_position*d_value, y_position*d_value, z_position*d_value,'.');
%% simulation
% LU decomp and bckward euler
tic;
% thermal simulation setup
t_stop = 12;
t_step = 1e-1;
t = 0:t_step:t_stop;
nstep = size(t,2);

% populate inputs along time axis
n_port = size(B,2);
n_np_source = n_port-n_power;
% Read watchh trace for power

% pw_trace contains 362 samples for 36 cores. First row is core number, 
% so we will skip it. Pick one core and pick n_step values from the trace
%watchh_trace=dlmread('pw_trace', ' ');
%pw_trace= watchh_trace(2:nstep+1,1)'; % 1 here is core number. Arbitrary
%u_vec = I'*pw_trace;
u_vec = I'*ones(1,nstep);
% set power to be 0 at time 0
u_vec(n_np_source+1:end,1) = zeros(n_power,1);
 save('GCB_matrix.mat','G','C','B','u_vec','t','t_step')
disp('Simulating, please wait, this may take several minutes...');
% thermal simulation without boundary node merge
[xres] = thermal_simulation_struct(G,C,B,u_vec,t_step);
%xres = thermal_static(G,C,B,u_vec);

xres = xres - 273.15;
sim_time = toc
%% plot
% power point plot
port_res = B'*xres;
%figure;
%plot(t,port_res(n_np_source+1,:));

% 3-D surface plot
n_fd = size(G,1);
xres_steady = xres(:,end);
xyz_data = cell(n_blk,1);
heat_map = zeros(1,1,1);
min_temp = 1000000;
max_temp = 0;
for i = 1:n_fd
    blk_idx = fd_to_blk{i}; % find fd_node belongs to which block
    xyz_data{blk_idx}(fd_to_xyz(i,1),fd_to_xyz(i,2),fd_to_xyz(i,3)) = xres_steady(i);
    %fd_to_xyz(i,:)
    heat_map(fd_to_xyz(i,1),fd_to_xyz(i,2),fd_to_xyz(i,3)) = xres_steady(i);
    if (xres_steady(i) > 0 && xres_steady(i) < min_temp)
        min_temp = xres_steady(i);
    end
    if (xres_steady(i) > 0 && xres_steady(i) > max_temp)
        max_temp = xres_steady(i);
    end
end

heat_map = 255 * (heat_map - min_temp) / (max_temp - min_temp);
for i = 1:n_blk
    xyz_data{i} = xyz_data{i}(blk_start{i}(1):end, blk_start{i}(2):end, blk_start{i}(3):end);
end
figure;
for i = 1:n_blk
    h = vol3d(d_value,dz_value, 'cdata',xyz_data{i},'texture','3D','start',blk_start{i});
    hold on;
end
hold off;
%h = vol3d(d_value,'cdata',xyz_data{2},'texture','3D','start',blk_start{2});
end
