function [k, K] = backward_pass(X, u, ts, N, l_cost_x, l_cost_xx, l_cost_u, l_cost_uu, l_cost_ux, f_dx, f_du, lamb)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[l_x, l_xx, l_u, l_uu, l_ux] = get_cost_derivative(X(:,2:end), u, ts, N, l_cost_x, l_cost_xx, l_cost_u, l_cost_uu, l_cost_ux);
[f_x, f_u] = get_system_derivative(X(:,2:end), u, ts, N, f_dx, f_du);

V_x = l_x(:,end);
V_xx = l_xx(:,:,end);

k = zeros(size(u,1), N);
K = zeros(size(u,1), size(X,1), N);

for i = (N-1):-1:1
   Q_x = l_x(:,i) + f_x(:,:,i)'*V_x;
   Q_u = l_u(:,i) + f_u(:,:,i)'*V_x;
   Q_xx = l_xx(:,:,i) + f_x(:,:,i)'*V_xx*f_x(:,:,i);
   Q_ux = l_ux(:,:,i) + f_u(:,:,i)'*V_xx*f_x(:,:,i);
   Q_uu = l_uu(:,:,i) + f_u(:,:,i)'*V_xx*f_u(:,:,i);
%    Q_uu_evalues = eig(Q_uu);
%    [Q_uu_evectors,Q_uu_evalues_aux] = eig(Q_uu);
   
   Q_uu(Q_uu<0)= 0;
   Q_uu = Q_uu + lamb;
   
%    diagonal = diag(1.0/Q_uu_evalues);
%    aux_1 = dot(diagonal,Q_uu_evectors');
%    Q_uu_inv = dot(Q_uu_evectors,aux_1);
   Q_uu_inv = inv(Q_uu);
   
   k(:,i) = - Q_uu_inv*Q_u;
   K(:,:,i) = -Q_uu_inv*Q_ux;
   
   V_x = Q_x -K(:,:,i)'*Q_uu*k(:,i);
   V_xx = Q_xx -K(:,:,i)'*Q_uu*K(:,:,i);
   
end
end

