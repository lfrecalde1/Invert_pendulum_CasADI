function [xp] = f_system(x, u, L)

%% DEFINITION SPACE DIMENTION
n_states = length(x);

%% CONSTANT VALUES
m0 = L(1);
m1 = L(2);
l1 = L(3);
I = L(4);
g = L(5);
b1 = L(6);
b2 = L(7);

%% UNMERGED VALUES
x1 = x(1);
th1 = x(2);
dx1 = x(3);
dth1 = x(4);

%% MATRICES OF THE SYSTEM

gama1 = I*m0 + I*m1 + l1^2*m1^2 + l1^2*m0*m1 - l1^2*m1^2*cos(th1)^2;

M_1 = [(m1*l1^2 + I)/gama1, -(l1*m1*cos(th1))/gama1;...
       -(l1*m1*cos(th1))/gama1, (m0 + m1)/gama1];
   
C = [0, -m1*l1*sin(th1)*dth1;...
     0,0];
 
G = [0;-m1*g*l1*sin(th1)];

T = [1;0];

F = [-b1, 0;...
      0, -b2];

%% REDUCE ORDER SYSTEM CREATION OF AXULIARY MATRICES
A = zeros(n_states,n_states);

A(1:2,3:4) = eye(2);

A(3:4,3:4) = -M_1*(C-F);

B = zeros(n_states,1);

B(3:4,1) = M_1*T;

E= zeros(n_states,1);

E(3:4,1) = -M_1*G;

xp = A*x+E+B*u;

end

