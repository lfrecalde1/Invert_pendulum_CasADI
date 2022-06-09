function [l_x, l_xx, l_u, l_uu, l_ux] = get_cost_derivative(x, u, ts, N, l_dx, l_dxx, l_du, l_duu, ldux)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
l_x = zeros(size(x,1),N);
l_xx = zeros(size(x,1),size(x,1),N);
l_u = zeros(size(u,1),N);
l_uu = zeros(size(u,1),size(u,1),N);
l_ux = zeros(size(u,1), size(x,1),N);

for k = 1:1:N
    l_x(:,k) = full(l_dx(x(:,k),u(:,k)))*ts;
    l_xx(:,:,k) = full(l_dxx(x(:,k),u(:,k)))*ts;
    l_u(:,k) = full(l_du(x(:,k),u(:,k)))*ts;
    l_uu(:,:,k) = full(l_duu(x(:,k),u(:,k)))*ts;
    l_ux(:,:,k) = full(ldux(x(:,k),u(:,k)))*ts;
end
end

