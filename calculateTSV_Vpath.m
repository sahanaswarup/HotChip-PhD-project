function [tkdd] = calculateTSV_Vpath(r_blks, fd_to_xyz, thermk_Bulk, TSV_locations_inCurrent_Blk, blks, fd_cur, fd_adj, d_value, dz_value, TSV, direction)
         
        [TSVIncludedCurrent, TSVsInGridCurrent, TSVxdimCurrent, TSVydimCurrent] = check_TSVs(fd_to_xyz(fd_cur,:), TSV_locations_inCurrent_Blk, TSV.radius, direction, d_value); 
 %       [TSVIncludedAdj, TSVsInGridAdj, TSVxdimAdj, TSVydimAdj] = check_TSVs(fd_to_xyz(fd_adj,:),TSV_locations_inCurrent_Blk, TSV.radius, inverse_direction(direction), d_value); 
       % if( TSVIncludedCurrent ~=0 && TSVIncludedAdj ~= 0)
       if( TSVIncludedCurrent ~=0 )

%             Rcu = (1/TSV.therm_k_Cu) * dz_value / ( pi* TSV.radius^2  * TSVxdimCurrent * TSVydimCurrent );
%             Rins = (1/TSV.therm_k_Ins) * dz_value / ( pi* TSV.radius^2 * ((1 + TSV.linear)^2 - 1) * TSVxdimCurrent * TSVydimCurrent );
%             RBulk = (1/thermk_Bulk) * dz_value / ( d_value * d_value - pi * TSV.radius^2 * (1 + TSV.linear)^2* TSVxdimCurrent * TSVydimCurrent );
            Rcu = (1/TSV.therm_k_Cu) * dz_value / ( pi* TSV.radius^2*(1-TSV.linear)^2  * TSVxdimCurrent * TSVydimCurrent );
            Rins = (1/TSV.therm_k_Ins) * dz_value / ( pi* TSV.radius^2 * (1-(1 - TSV.linear)^2) * TSVxdimCurrent * TSVydimCurrent );
            RBulk = (1/thermk_Bulk) * dz_value / ( d_value * d_value - pi * TSV.radius^2 * TSVxdimCurrent * TSVydimCurrent );

            R_vertical = 1/(1/Rcu + 1/Rins + 1/RBulk);
            tkdd = 1/(R_vertical * d_value^2 * dz_value); 
%               therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
%               tkdd = therm_k_eq/dz_value/dz_value/2; 
%         elseif ( TSVIncludedCurrent ~=0 && TSVIncludedAdj == 0)
%             %display('calculateTSV_Vpath2');
%             Rcu = (1/TSV.therm_k_Cu) * d_value / ( pi* TSV.radius^2 * (1 - TSV.linear^2) * TSVxdimCurrent * TSVydimCurrent ) /2 ;
%             Rins = (1/TSV.therm_k_Ins) * d_value / ( pi* TSV.radius^2 * TSV.linear^2 * TSVxdimCurrent * TSVydimCurrent ) / 2;
%             RBulk = blks{r_blks(1)}.therm_k * d_value / ( d_value * d_value - pi * TSV.radius^2 * TSVxdimCurrent * TSVydimCurrent ) / 2;
%             R_vertical_currenthalf = 1/(1/Rcu + 1/Rins + 1/RBulk);
%             R_vertical_adjacenthalf = blks{r_blks(2)}.therm_k * d_value /(d_value * d_value);
%             R_vertical = R_vertical_currenthalf + R_vertical_adjacenthalf; 
%             tkdd = R_vertical/d_value^3; 
        else  %bulk material only in vertical direction
             %therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
             % tkdd = therm_k_eq/dz_value/dz_value/2
              therm_k_eq = min([blks{r_blks(1)}.therm_k, blks{r_blks(2)}.therm_k]);
              tkdd = therm_k_eq/dz_value/dz_value;
        %      tkdd = (blks{r_blks(1)}.therm_k/dz_value/(dz_value/2) * blks{r_blks(2)}.therm_k/dz_value/(dz_value/2))/(blks{r_blks(1)}.therm_k/dz_value/(dz_value/2) + blks{r_blks(2)}.therm_k/dz_value/(dz_value/2)); 
        end 
