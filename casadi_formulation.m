function [f,fd,f_x,f_u, l, l_u, l_uu, l_x, l_xx, l_ux, l_f, l_fx, l_fxx] = casadi_formulation(L, ts)

addpath('/home/fer/casadi-linux-matlabR2014b-v3.4.5');
import casadi.*;

%% CONSTANT VALUES
m0 = L(1);
m1 = L(2);
l1 = L(3);
I = L(4);
g = L(5);
b1 = L(6);
b2 = L(7);

%% DEFINE SYSTEM VARIABLES
u1 = SX.sym('u1');
controls = [u1];
n_controls = length(controls);

x1 = SX.sym('x1'); 
th1 = SX.sym('th1'); 
dx1 = SX.sym('dx1'); 
dth1 = SX.sym('dth1'); 
states = [x1;th1;dx1;dth1];

x1_d = SX.sym('x1_d'); 
th1_d = SX.sym('th1_d'); 
dx1_d = SX.sym('dx1_d'); 
dth1_d = SX.sym('dth1_d'); 
states_d = [x1_d;th1_d;dx1_d;dth1_d];

n_states = length(states);
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
A = [zeros(2,2),eye(2);...
     zeros(2,2),-M_1*(C-F)];

B = [zeros(2,1);...
    M_1*T];

E = [zeros(2,1);...
    -M_1*G];

%% SYSTEM DEFINITION
rhs = A*states+E+B*controls;
f = Function('f',{states,controls},{rhs});
% compute solution symbolically


%% RUNGE KUTA 4
k1 = f(states, controls);   % new
k2 = f(states + ts/2*k1, controls); % new
k3 = f(states + ts/2*k2, controls); % new
k4 = f(states + ts*k3, controls); % new
x_k = states +ts/6*(k1 +2*k2 +2*k3 +k4); % new
fd=Function('fd',{states,controls},{x_k});


%% PARTIAL DERIVATIES SYSTEM
f_dx = jacobian(rhs,states);
f_du = jacobian(rhs,controls);
f_x = Function('f_x',{states,controls},{f_dx});
f_u = Function('f_u',{states,controls},{f_du});

%% COST FUNCTION
error =states_d-states;

Q = [50,0,0,0;...
     0,50,0,0;...
     0,0,0,0;...
     0,0,0,0];
 
R =1*eye(1);
l_cost = (1/2)*(states'*Q*states)+(1/2)*(controls'*R*controls);
l_final = (1/2)*(states'*Q*states);

l = Function('l',{states, controls},{l_cost});
l_f = Function('l',{states},{l_final});

l_u = Function('l_u',{states, controls},{jacobian(l_cost,controls)});
l_uu = Function('l_uu',{states, controls},{hessian(l_cost,controls)});

l_u_aux = gradient(l_cost,controls);
l_ux = Function('l_ux',{states, controls},{jacobian(l_u_aux, states)});

l_x = Function('l_x',{states, controls},{jacobian(l_cost,states)});
l_xx = Function('l_xx',{states, controls},{hessian(l_cost,states)});

l_fx = Function('l_fx',{states},{jacobian(l_final,states)});
l_fxx = Function('l_fxx',{states},{hessian(l_final,states)});
end
