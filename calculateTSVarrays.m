function [tkdd] = calculateTSVarrays(r_blks, fd_to_xyz, thermk_Bulk, TSV_locations_inCurrent_Blk, blks, fd_cur,fd_adj, d_value, dz_value, TSV, direction)
         
%         beta1 = 0.027943775072880; 
%         beta2 =  0.956434460745102; 
          beta1 =   -0.007264868908017;
          beta2 = 0.086062726337019;
          beta3 = 0.869256033848894; 
        [TSVIncludedCurrent, TSVsInGridCurrent, TSVxdimCurrent, TSVydimCurrent] = check_TSVs(fd_to_xyz(fd_cur,:), TSV_locations_inCurrent_Blk, TSV.radius, direction, d_value);
        [TSVIncludedAdj, TSVsInGridAdj, TSVxdimAdj, TSVydimAdj] = check_TSVs(fd_to_xyz(fd_adj,:),TSV_locations_inCurrent_Blk, TSV.radius, inverse_direction(direction), d_value); 
        TSVarraydimX = TSVxdimCurrent + TSVxdimAdj;
        TSVarraydimY = TSVydimCurrent + TSVydimAdj;
        
        %TSV_space is a fixed variable
%         if(TSVIncludedCurrent ~= 0) %check if current grid has TSVs
%             TSVsInGrid = TSVsInGridCurrent;
%         elseif(TSVIncludedAdj ~= 0) %check if adjacent grid has TSVs
%             TSVsInGrid = TSVsInGridAdj;
%         else     %no TSVs in both current grid and adjacent grid
%             TSVsInGrid = []; 
%         end
        TSVsInGrid = [TSVsInGridCurrent; TSVsInGridAdj];
        if (size(TSVsInGrid,1) <= 1)
            TSV_space = d_value - 2*TSV.radius*(1+TSV.linear);
            TSVarray_parallel = 1;
            TSVarray_vertical = 1;
        else
            switch direction
                case {'left','right'}
                    if(abs(TSVsInGrid(1,1) - TSVsInGrid(2,1)) ~= 0)
                        TSVcenterTocenter = abs(TSVsInGrid(1,1) - TSVsInGrid(2,1));
                    elseif (abs(TSVsInGrid(1,2) - TSVsInGrid(2,2)) ~= 0)
                        TSVcenterTocenter = abs(TSVsInGrid(1,2) - TSVsInGrid(2,2));
                    else
                        TSVcenterTocenter = d_value; 
                    end
                   TSV_space = TSVcenterTocenter - 2*TSV.radius*(1+TSV.linear);                     
                   %TSV_space = abs (TSVsInGrid(1,2) - TSVsInGrid(2,2)) - 2*TSV.radius*(1+TSV.linear); 
                   TSVarray_parallel = TSVarraydimY;
                   TSVarray_vertical = TSVarraydimX; 
                case {'front','back'}
                    if(abs(TSVsInGrid(1,1) - TSVsInGrid(2,1)) ~= 0)
                        TSVcenterTocenter = abs(TSVsInGrid(1,1) - TSVsInGrid(2,1));
                    elseif (abs(TSVsInGrid(1,2) - TSVsInGrid(2,2)) ~= 0)
                        TSVcenterTocenter = abs(TSVsInGrid(1,2) - TSVsInGrid(2,2));
                    else
                        TSVcenterTocenter = d_value; 
                    end
                   TSV_space = TSVcenterTocenter - 2*TSV.radius*(1+TSV.linear); 
                   %TSV_space = abs(TSVsInGrid(1,1) - TSVsInGrid(2,1)) - 2*TSV.radius*(1+TSV.linear); 
                   TSVarray_parallel = TSVarraydimX;
                   TSVarray_vertical = TSVarraydimY; 
            end
        end
        if( TSVIncludedCurrent ~=0 && TSVIncludedAdj ~=0)
            TSVarray_parallel_dim = TSVarray_parallel/2;
        else
            TSVarray_parallel_dim = TSVarray_parallel;
        end

        %thermk_Bulk is the thermal conductivty of the bulk silicon for
        %TSVs by its definition
        if( TSVIncludedCurrent ~=0 || TSVIncludedAdj ~=0)
%       calculate the lateral thermal resistance of TSV array
%       C(i,i) = (TSV.rho*TSV.cp);
             %display('fuck'); 
%             theta = beta1 * TSV_space / (TSV.radius) + beta2; 
             theta = beta1 * (TSV_space / (TSV.radius))^2 + beta2 * (TSV_space / (TSV.radius)) + beta3;
             RTSV = -2*(1/TSV.therm_k_Ins)*log(1-TSV.linear)/(pi*dz_value) + (1/TSV.therm_k_Cu)*pi/(2*dz_value);
             %Rspace = (1/thermk_Bulk)* 2*TSV.radius*(1+TSV.linear) /(TSV_space *d_value); 
             Rspace = (1/thermk_Bulk)* 2*TSV.radius /(TSV_space *dz_value); 
             RTSVarray = RTSV * Rspace /(RTSV + Rspace) / TSVarray_parallel_dim * TSVarray_vertical *  theta  ;
             RSi_Serial = (1/thermk_Bulk)*(d_value - 2*TSV.radius * TSVarray_vertical )/(d_value*dz_value);
             tkdd = 1/((RSi_Serial + RTSVarray  )*(d_value^2*dz_value)); 
%%%%%%%%%%%%%%%%%%%Xu's method%%%%%%%%%%%%%%%%%5
%             kl = 1/0.4465;
%             RTSV = -(1/TSV.therm_k_Ins)*log(1-TSV.linear)/(2*pi*kl*dz_value);
%             RTSVarray = RTSV/TSVarray_parallel_dim * TSVarray_vertical;
%             RSi_Serial = (1/thermk_Bulk)*(d_value - 2*TSV.radius * TSVarray_vertical )/(d_value*dz_value);
%             tkdd = 1/((RSi_Serial + RTSVarray  )*(d_value^2*dz_value)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
%                 tkdd = therm_k_eq/d_value/d_value/2;
%                  therm_k_eq = blks{r_blks(1)}.therm_k * (d_value^2 - pi * TSV.radius^2*(1+TSV.linear)^2 * TSVarray_parallel_dim * TSVarray_vertical )/d_value^2 + ...
%                    TSV.therm_k_Cu * (pi * TSV.radius^2 *TSVarray_parallel_dim * TSVarray_vertical )/d_value^2+...
%                    TSV.therm_k_Ins * (pi * TSV.radius^2*((1+TSV.linear)^2-1) * TSVarray_parallel_dim * TSVarray_vertical )/d_value^2;
%                tkdd = therm_k_eq/d_value^2; 
        else  %this grid and the one adjacent to it have no TSVs inside
              therm_k_eq =  blks{r_blks(1)}.therm_k+blks{r_blks(2)}.therm_k;
              tkdd = therm_k_eq/d_value/d_value/2;                 
        end 
