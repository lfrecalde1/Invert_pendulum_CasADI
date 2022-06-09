function [x_new,u_new] = forward_pass(x, u, ts, N, f, k, K)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x_new = zeros(size(x,1),N+1);
x_new(:,1) = x(:,1);
u_new = zeros(size(u,1), N);
for i = 1:1:N
   u_new(:,i) = u(:,i) + k(:,i) + K(:,:,i)*(x_new(:,i)-x(:,i));
   x_new(:,i+1) = forward_simulate( x_new(:,i), u_new(:,i), ts, f);
end
end

