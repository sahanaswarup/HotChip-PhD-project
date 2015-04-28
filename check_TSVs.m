function [IsIncluded, TSV_locsInGrid, TSV_xdim, TSV_ydim]  = check_TSVs(xyz, TSV_locs, TSV_radius, direction, d_value)

 if isempty(TSV_locs)
     IsIncluded = 0;
     TSV_xdim = 0; 
     TSV_ydim = 0; 
     TSV_locsInGrid = []; 
 else
     TSV_xdim = 0; 
     TSV_ydim = 0; 
     IsIncluded = 0;
     TSV_locsInGrid = []; 
     %find out which of the TSVs are in the grid centered at xyz
     % 'left' indicates that only the TSVs in the left side of the cell
     % need to be considered. 
     Econst = 10000;
     left_limit   = xyz(1) - 1;
     right_limit  = xyz(1);
     back_limit   = xyz(2) + 1;
     front_limit  = xyz(2);
     center_x     = xyz(1) - 1/2;
     center_y     = xyz(2) + 1/2;
     switch direction
         case 'left'
             for k=1:length(TSV_locs)
                  if((TSV_locs{k}(1) - TSV_radius)/d_value >= left_limit && (TSV_locs{k}(1) + TSV_radius)/d_value <= center_x && ...
                     (TSV_locs{k}(2) - TSV_radius)/d_value >= front_limit && (TSV_locs{k}(2) + TSV_radius)/d_value <= back_limit ...
                     || ( ((TSV_locs{k}(1)/d_value) - center_x)*Econst  == 0 && ...
                          (TSV_locs{k}(2) + TSV_radius)/d_value <= back_limit && ...
                          (TSV_locs{k}(2) - TSV_radius)/d_value >= front_limit))         
                      IsIncluded = 1;
                      TSV_locsInGrid  = [TSV_locsInGrid; TSV_locs{k}(1),TSV_locs{k}(2)];
                  end
             end
         case 'right'
            for k=1:length(TSV_locs)
                  if ((TSV_locs{k}(1) - TSV_radius)/d_value >= center_x && (TSV_locs{k}(1) + TSV_radius)/d_value <= right_limit && ...
                      (TSV_locs{k}(2) - TSV_radius)/d_value >= front_limit && (TSV_locs{k}(2) + TSV_radius)/d_value <= back_limit ...
                       || (((TSV_locs{k}(1)/d_value)- center_x)*Econst  == 0 && ...
                           (TSV_locs{k}(2) + TSV_radius)/d_value <= back_limit && ...
                           (TSV_locs{k}(2) - TSV_radius)/d_value >= front_limit))         
                      IsIncluded = 1;
                      TSV_locsInGrid  = [TSV_locsInGrid; TSV_locs{k}(1),TSV_locs{k}(2)];
                  end
             end
         case 'front'
             for k=1:length(TSV_locs)
                  if((TSV_locs{k}(1) - TSV_radius)/d_value >= left_limit && (TSV_locs{k}(1) + TSV_radius)/d_value <= right_limit && ...
                     (TSV_locs{k}(2) - TSV_radius)/d_value >= front_limit && (TSV_locs{k}(2) + TSV_radius)/d_value <= center_y...
                      || (((TSV_locs{k}(2)/d_value)- center_y)*Econst  == 0 && ...
                           (TSV_locs{k}(1) + TSV_radius)/d_value <= right_limit && ...
                           (TSV_locs{k}(1) - TSV_radius)/d_value >= left_limit))         
                      IsIncluded = 1;
                      TSV_locsInGrid  = [TSV_locsInGrid; TSV_locs{k}(1),TSV_locs{k}(2)];
                  end
             end
         case 'back'
             for k=1:length(TSV_locs)
                  if((TSV_locs{k}(1) - TSV_radius)/d_value >= left_limit && (TSV_locs{k}(1) + TSV_radius)/d_value <= right_limit && ...
                     (TSV_locs{k}(2) - TSV_radius)/d_value >= center_y && (TSV_locs{k}(2) + TSV_radius)/d_value <= back_limit ...
                     || (((TSV_locs{k}(2)/d_value)- center_y)*Econst  == 0 && ...
                          (TSV_locs{k}(1) + TSV_radius)/d_value <= right_limit && ...
                          (TSV_locs{k}(1) - TSV_radius)/d_value >= left_limit))         
                      IsIncluded = 1;
                      TSV_locsInGrid  = [TSV_locsInGrid; TSV_locs{k}(1),TSV_locs{k}(2)];
                  end
             end
         case {'up','down','all'}
             for k=1:length(TSV_locs)
                  if ((TSV_locs{k}(1) - TSV_radius)/d_value >= left_limit && (TSV_locs{k}(1) + TSV_radius)/d_value <= right_limit && ...
                      (TSV_locs{k}(2) - TSV_radius)/d_value >= front_limit && (TSV_locs{k}(2) + TSV_radius)/d_value <= back_limit ...
                     )         
                      IsIncluded = 1;
                      TSV_locsInGrid  = [TSV_locsInGrid; TSV_locs{k}(1),TSV_locs{k}(2)];
                  end
             end
     end
     
     if isempty(TSV_locsInGrid)  %when no TSV is found in the grid
         return; 
     end
     %find out the dimension of TSV
     prev_xlocs = -1; 
     sortedxlocs = sort(TSV_locsInGrid(:,1));
     prev_ylocs = -1; 
     sortedylocs = sort(TSV_locsInGrid(:,2));
     
     for m=1:size(TSV_locsInGrid, 1)
         xlocs =  sortedxlocs(m); 
         if(xlocs ~= prev_xlocs )
             xlocs = sortedxlocs(m);
             TSV_xdim = TSV_xdim + 1; 
         end
         prev_xlocs = xlocs; 
     end            
     for m=1:size(TSV_locsInGrid, 1)
         ylocs = sortedylocs(m); 
         if(ylocs ~= prev_ylocs)
             ylocs = sortedylocs(m);
             TSV_ydim = TSV_ydim + 1; 
         end
         prev_ylocs = ylocs; 
     end 
     


 end