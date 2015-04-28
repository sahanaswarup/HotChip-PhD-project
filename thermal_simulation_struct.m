function [xres] = thermal_simulation_struct(G,C,B,u_vec,h)

n_t_step = size(u_vec,2); % number of timesteps
N = n_t_step-1;

% DC
wn1 = B*u_vec;
xn = G\wn1(:,1);
xdc = xn;

% backward
% trans
left = G + 1/h*C;
tic
[L,U,P,Q] = lu(left);
toc
right = 1/h*C;
for i = 1:N
    b = right*xn + wn1(:,i+1);
    b=sparse(b);
    xn1 = Q*(U\(L\(P*b)));  
    xres(:,i) = xn1;
    xn = xn1;
end;

xres = [xdc xres];

