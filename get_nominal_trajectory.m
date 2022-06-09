function [x] = get_nominal_trajectory(x0, u, ts, N, f)
x = zeros(size(x0,1),N+1);
x(:,1) = x0;

for k = 1:1:N
   x(:,k+1) = forward_simulate( x(:,k), u(:,k), ts, f);
end

end

