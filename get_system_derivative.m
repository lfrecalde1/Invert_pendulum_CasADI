function [f_x, f_u] = get_system_derivative(x, u, ts, N, f_dx, f_du)
i = size(x,1);
j = size(u,1);
f_x = zeros(i, i, N);
f_u = zeros(i, j, N);

I = eye(i);
for k =1:N
   A = full(f_dx(x(:,k),u(:,k)));
   B = full(f_du(x(:,k),u(:,k)));
   f_x(:,:,k) = I + A*ts;
   f_u(:,:,k) = B*ts;
end
end

