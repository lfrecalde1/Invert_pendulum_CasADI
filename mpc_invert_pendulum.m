function [f,solver,args] = mpc_invert_pendulum(bounded, N, L, ts)

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

%% Limita variables
ul_max = bounded(1);
ul_min = bounded(2);

%% DEFINE SYSTEM VARIABLES
u1 = SX.sym('u1');
controls = [u1];
n_controls = length(controls);

x1 = SX.sym('x1'); 
th1 = SX.sym('th1'); 
dx1 = SX.sym('dx1'); 
dth1 = SX.sym('dth1'); 
states = [x1;th1;dx1;dth1];


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

U = SX.sym('U',n_controls,N);
P = SX.sym('P',n_states + N*(n_states));
%% vector que incluye el vector de estados y la referencia
X = SX.sym('X',n_states,(N+1));

%% Vector que representa el problema de optimizacion
g = [];  % restricciones de estados del problema  de optimizacion

%%EMPY VECTOR ERRORS
he = [];

%% EMPY VECTOR CONTROL VALUES
u = [];
%% INITIAL CONDITION ERROR
g = [g;X(:,1)-P(1:4)]; % initial condition constraints

%% Definicon del bucle para optimizacion a los largo del tiempo
for k = 1:N
    st = X(:,k);  con = U(:,k);

    %% Funcion costo a minimizar 
    he = [he;X(1:2,k)-P(4*k+1:4*k+2)];
    u = [u;con];
    %obj = obj+(st-P(4*k+1:4*k+4))'*Q*(st-P(4*k+1:4*k+4)) + con'*R*con;
    
    %% Actualizacion del sistema usando Euler runge kutta
    st_next = X(:,k+1);
    k1 = f(st, con);   % new 
    k2 = f(st + ts/2*k1, con); % new
    k3 = f(st + ts/2*k2, con); % new
    k4 = f(st + ts*k3, con); % new
    st_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4); % new 
    
    %% Restricciones del sistema se =basan en el modelo del sistema
    g = [g;st_next-st_next_RK4]; 
end

%% Cost final 
Q = 50*eye(size(he,1));
R = 1*eye(size(u,1));

%% FINAL COST
obj = he'*Q*he+u'*R*u;
% se crea el vector de desiscion solo de una columna
OPT_variables = [reshape(X,4*(N+1),1);reshape(U,1*N,1)];

nlprob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 1;
opts.ipopt.acceptable_tol =1e-4;
opts.ipopt.acceptable_obj_change_tol = 1e-4;

solver = nlpsol('solver', 'ipopt', nlprob,opts);

args = struct;

args.lbg(1:4*(N+1)) = 0;  %-1e-20  %Equality constraints
args.ubg(1:4*(N+1)) = 0;  %1e-20   %Equality constraints

args.lbx(1:4:4*(N+1),1) = -inf; %state x lower bound
args.ubx(1:4:4*(N+1),1) = inf;  %state x upper bound

args.lbx(2:4:4*(N+1),1) = -inf; %state y lower bound
args.ubx(2:4:4*(N+1),1) = inf;  %state y upper bound

args.lbx(3:4:4*(N+1),1) = -inf; %state z lower bound
args.ubx(3:4:4*(N+1),1) = inf;  %state z upper bound

args.lbx(4:4:4*(N+1),1) = -inf; %state theta lower bound
args.ubx(4:4:4*(N+1),1) = inf;  %state theta upper bound

args.lbx(4*(N+1)+1:1:4*(N+1)+1*N,1) = ul_min;  %
args.ubx(4*(N+1)+1:1:4*(N+1)+1*N,1) = ul_max;  %


end