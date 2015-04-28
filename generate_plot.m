% generating file that Hai's code takes as input

function  generate_plot(x_y_dim, z_dim, num_layers, d_value, dz_value, ...
                        gap_bw_layers, tsv_option, tsv_file, tsv_r, tsv_t, ...
                        filename)

    scale         = 1;
    x_start       = 0.25 * x_y_dim / scale;
    y_start       = 0.25 * x_y_dim / scale;
    x_dim         = x_y_dim / scale;
    y_dim         = x_y_dim / scale;
    
    x_start_hs    = 0;
    y_start_hs    = 0;
    x_dim_hs      = (1.5 * x_y_dim) / scale;
    y_dim_hs      = (1.5 * x_y_dim) / scale;
    
    z_dim         = z_dim / scale;
    gap_bw_layers = gap_bw_layers / scale;
    tsv_r         = tsv_r / scale;
    tsv_t         = tsv_t / scale;
    d_value       = d_value / scale;
    dz_value      = dz_value / scale;
    num_cells     = x_dim / d_value;
    
    fp = fopen(filename, 'w');

    fprintf(fp, 'scale 10000\n');
    fprintf(fp, 'd_value %.2f\n', d_value);
    fprintf(fp, 'dz_value %.2f\n', dz_value);
    fprintf(fp, 't_amb 20\n');
    fprintf(fp, 'TSV_therm_k_Cu  400\n');
    fprintf(fp, 'TSV_rho_Cu 8700\n');
    fprintf(fp, 'TSV_cp_Cu 385\n');
    fprintf(fp, 'TSV_therm_k_Ins  1.38\n');
    fprintf(fp, 'TSV_rho_Ins 2650\n');
    fprintf(fp, 'TSV_cp_Ins 100\n');
    fprintf(fp, 'TSV_radius %.3f\n', tsv_r);
    fprintf(fp, 'TSV_linear %.3f\n\n\n', tsv_t/tsv_r);

    z_start = 0;
    for k=1:num_layers
        fprintf(fp, 'block_begin\n');
        fprintf(fp, 'coordinate %.8f %.8f %.8f\n', x_start, y_start, z_start);
        fprintf(fp, 'dim %.3f %.3f %.3f\n', x_dim, y_dim, z_dim/2);
        fprintf(fp, 'therm_k 160\n');
        fprintf(fp, 'he 0 0 0 0 0 0\n');
        fprintf(fp, 'rho 2300\n');
        fprintf(fp, 'cp 700\n');
        fprintf(fp, 'channel 0\n');
        if tsv_option == 1
            for i=1:num_cells
                for j=1:num_cells
                    x_pos = x_start + (i - 1/2) * d_value;
                    y_pos = y_start + (j - 1/2) * d_value;
                    fprintf(fp, 'TSV_locations %.8f %.8f\n', x_pos, y_pos);
                end
            end
        else
            tsv_fp = csvread(tsv_file);
            for i = 1:size(tsv_fp, 1)
                cell_x_id = tsv_fp(i,1);
                cell_y_id = tsv_fp(i,2);
                if (cell_x_id < 1 || cell_x_id > num_cells || ...
                    cell_y_id < 1 || cell_y_id > num_cells)
                    continue
                end
                x_pos = x_start + (cell_x_id - 1/2) * d_value;
                y_pos = y_start + (cell_y_id - 1/2) * d_value;
                layer = tsv_fp(i,3);
                if layer == k
                    fprintf(fp, 'TSV_locations %.8f %.8f\n', x_pos, y_pos);
                end
            end
        end
        fprintf(fp, 'block_end\n\n\n');
        % metal layer
        z_start = z_start + z_dim/2;
        fprintf(fp, 'block_begin\n');
        fprintf(fp, 'coordinate %.8f %.8f %.8f\n', x_start, y_start, z_start);
        fprintf(fp, 'dim %.3f %.3f %.3f\n', x_dim, y_dim, z_dim/2);
        fprintf(fp, 'therm_k 160\n');
        fprintf(fp, 'he 0 0 0 0 0 0\n');
        fprintf(fp, 'rho 2300\n');
        fprintf(fp, 'cp 700\n');
        fprintf(fp, 'channel 0\n');
        if tsv_option == 1
            for i=1:num_cells
                for j=1:num_cells
                    x_pos = x_start + (i - 1/2) * d_value;
                    y_pos = y_start + (j - 1/2) * d_value;
                    fprintf(fp, 'TSV_locations %.8f %.8f\n', x_pos, y_pos);
                end
            end
        else
            tsv_fp = csvread(tsv_file);
            for i = 1:size(tsv_fp, 1)
                cell_x_id = tsv_fp(i,1);
                cell_y_id = tsv_fp(i,2);
                if (cell_x_id < 1 || cell_x_id > num_cells || ...
                    cell_y_id < 1 || cell_y_id > num_cells)
                    continue
                end
                x_pos = x_start + (cell_x_id - 1/2) * d_value;
                y_pos = y_start + (cell_y_id - 1/2) * d_value;
                layer = tsv_fp(i,3);
                if layer == k
                    fprintf(fp, 'TSV_locations %.8f %.8f\n', x_pos, y_pos);
                end
            end
        end
        fprintf(fp, 'block_end\n\n\n');
        z_start = z_start + z_dim/2 + gap_bw_layers;
    end
    
    fprintf(fp, 'block_begin\n');
    fprintf(fp, 'coordinate %.3f %.3f %.3f\n', x_start_hs, y_start_hs, z_start);
    fprintf(fp, 'dim %.3f %.3f %.3f\n', x_dim_hs, y_dim_hs, 4 * dz_value); % This z_dim should be thinkness of heat sink which can be exposed in UI
    fprintf(fp, 'therm_k 160\n');
    fprintf(fp, 'he 0 6000.000000 0 0 0 0\n');
    fprintf(fp, 'rho 2700\n');
    fprintf(fp, 'cp 900\n');
    fprintf(fp, 'channel 0\n');
    fprintf(fp, 'block_end\n');
    fclose(fp);
end

