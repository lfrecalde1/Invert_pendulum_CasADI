%% Program to verifiy parcial derivaties using Casadi
clc, clear all, close all;

%% Path Casadi Package
addpath('/home/fer/casadi-linux-matlabR2014b-v3.4.5');
import casadi.*;


%% Define control actiones
u1 = SX.sym('u1'); 
u2 = SX.sym('u2');
U = [u1;u2];

%% Define cost over time
l_cost = (1/2)*(U'*U);
l = Function('l',{U},{l_cost});
l_u = Function('l_u',{U},{jacobian(l_cost,U)});
l_uu = Function('l_uu',{U},{hessian(l_cost,U)});

full(l_uu([10,10]))