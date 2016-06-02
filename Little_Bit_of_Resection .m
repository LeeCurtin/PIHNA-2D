% y = 75;
%       y_count = y;
%       y_max = 147;
%       x = 130;
%       x_count = x;
%       x_max = 185;
%       ind_res = [];
%       for i = sub2ind([y_max,x_max],y,x):N
%           if i >= sub2ind([y_max,x_max],y,x_count) && i <=sub2ind([y_max,x_max],y_max,x_count)
%             c(i) = 0;
%             h(i) = 0;
%             v(i) = 0;
%             a(i) = 0;
%             n(i) = 0;
%             ind_res = [ind_res;i];
%           end
%           
%           if mod(i,y_max) == 0
%             x_count = x_count+1;
%           end
%       end
%       
%       %ISCHEMIC EVENT
%       x_isch = 110;
%       y_isch = 75;
%       x_count = x_isch;
%       for i = sub2ind([y_max,x_max],y_isch,x_isch):N
%           if i >= sub2ind([y_max,x_max],y_isch,x_count) && i <=sub2ind([y_max,x_max],y_max,x_count)
%             v(i) = 0;
%           end
%           
%           if mod(i,y_max) == 0
%             x_count = x_count+1;
%           end
%       end
%       
%       
%       [Dweight1, Dweight2, Origin1, Origin2, Terminus1, Terminus2, ~, ~, Slice] ...
%           = set_resection_interfaces(Slice, ind_res);
%       
%       for i = 1:ind_nz
%           I = ind_nz(i);
%           if  ~ismember(I-1,ind_nz) && ~ismember(I+1,ind_nz) && ~ismember(I-n1,ind_nz) && ~ismember(I+n1,ind_nz)
%             c(i) = 0;
%             h(i) = 0;
%             v(i) = 0;
%             a(i) = 0;
%             n(i) = 0; 
%             Slice(i) = 0; %Deletes isolated tissue cells
%           end
%       end
%       
%       for i = 1:3
%           quick_fix = sub2ind([y_max,x_max],121+i,129);
% 
%           c(quick_fix) = 0;
%           h(quick_fix) = 0;
%           v(quick_fix) = 0;
%           a(quick_fix) = 0;
%           n(quick_fix) = 0;
%           Slice(i) = 0;
%       end
%       
%       Resection = 1; %Resection is over
%        
%     end