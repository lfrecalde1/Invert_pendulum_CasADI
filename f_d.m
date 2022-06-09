function [x] = f_d(x,u,L,ts)
k1 = f_system(x, u, L);   % new
k2 = f_system(x + ts/2*k1, u, L); % new
k3 = f_system(x + ts/2*k2, u, L); % new
k4 = f_system(x + ts*k3, u, L); % new
x = x +ts/6*(k1 +2*k2 +2*k3 +k4); % new
end

