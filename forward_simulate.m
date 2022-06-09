function [x] = forward_simulate(x, u, ts, f)
%UNTITLED2 Summary of this function goes here
x = x + full(f(x,u))*ts;
end

