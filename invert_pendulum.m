%% PROGRAM INVERT PENDULUM OVER CAR
%% CLEAN VARIABLES
clc, clear all, close all;

%% SET TIME VARIABLES
to = 0;
t_final = 30;
ts = 0.01;
t = (to:ts:t_final);
T =10;
t_aux = (to:ts:T);
%% SET SYSTEM CONSTANTS
g = 9.8;
m0 = 0.48;
m1 = 0.16;
L1 = 1;
l1 = L1/2;
I = (m1*l1^2)/12;
b1 = 1;
b2 = 0.0218;
%% GENERAL VECTOR OF CONSTAN VALUES
L = [m0; m1; l1; I; g; b1; b2];

%% CREATION EMPTY MATRIX FOR THE EVOLUTION OF THE SYSTEM
x = zeros(4,length(t)+1);
x(1,1) = 5; %% INITIAL POSITION OF THE CAR
x(2,1) = 180*(pi/180); %% INITIAL ANGULAR DISPLACEMENT OF THE PENDULUM
x(3,1) = 0; %% INITIAL VELOCITY CAR
x(4,1) = 0; %% INITIAL ANGULAR VELOCITUY OF THE PENDULUM

%% INITIALIZATION CONTROL VECTOR
u = ones(1,length(t));
U = u;
[f_c, f_discrete, f_dx, f_du, l_cost, l_cost_u, l_cost_uu, l_cost_x, l_cost_xx, l_cost_ux, l_final, l_final_x, l_final_xx] = casadi_formulation(L, ts);

%% GAIN MATRICES
Q=diag([1 1 1 1]);
R=1;

%% DESIRED SIGNAL
x1d = 0*sin(0.2*t);
x2d = 0*ones(1,length(t));
x3d = 0*ones(1,length(t));
x4d = 0*ones(1,length(t));

xd = [x1d;...
      x2d;...
      x3d;...
      x4d];
  
N = length(t_aux);
for i=1:length(t) 
    A(:,:,i) = full(f_dx(x(:,i),u(:,i)));
    B(:,:,i) = full(f_du(x(:,i),u(:,i)));
    x(:,i+1) = x(:,i) + full(f_c(x(:,i),u(:,i)))*ts;
end

X = get_nominal_trajectory(x(:,1), U, ts, N, f_c);
J_old = realmax;
lamb = 1;
lamb_max =10;
lamb_factor = 2;
tol = 0.01;
for itr=1:1:100
    [k, K] = backward_pass(X, U, ts, N, l_cost_x, l_cost_xx, l_cost_u, l_cost_uu, l_cost_ux, f_dx, f_du, lamb);
    [x_new,u_new] = forward_pass(X, U, ts, N, f_c, k, K);
    J_new = get_total_cost(X, U, ts, N, l_cost, l_final);
    if J_new < J_old
        X = x_new;
        U = u_new;
        lamb = (lamb)/lamb_factor;
        disp("Disminuye");
        if (abs(J_old-J_new)<tol)
            disp("Limit tolerance");
            break
        end
    else
        lamb = lamb*lamb_factor;
        disp("No disminuye")
        if lamb > lamb_max
            break
        end
    end
    J_old =J_new;
    J_new
end

figure(1)
for k=1:10:N
    drawpend(X(:,k),m1,m0,L1);
    pause(ts)
end
%% Parameters fancy plots
% define plot properties
lw = 2; % linewidth 1
lwV = 2; % linewidth 2
fontsizeLabel = 9; %11
fontsizeLegend = 9;
fontsizeTicks = 9;
fontsizeTitel = 9;
sizeX = 900; % size figure
sizeY = 300; % size figure

% color propreties
C1 = [246 170 141]/255;
C2 = [51 187 238]/255;
C3 = [0 153 136]/255;
C4 = [238 119 51]/255;
C5 = [204 51 17]/255;
C6 = [238 51 119]/255;
C7 = [187 187 187]/255;
C8 = [80 80 80]/255;
C9 = [140 140 140]/255;
C10 = [0 128 255]/255;
C11 = [234 52 89]/255;
C12 = [39 124 252]/255;
C13 = [40 122 125]/255;
C14 = [252 94 158]/255;
C15 = [244 171 39]/255;
C16 = [100 121 162]/255;
C17 = [255 0 0]/255;

C18 = [242 144 182]/255;
C19 = [229 23 101]/255;
C20 = [129 231 174]/255;
C21 = [42 142 86]/255;
C22 = [139 177 221]/255;
C23 = [64 109 159]/255;
C24 = [233 230 129]/255;
C25 = [178 173 12]/255;

figure('Position', [15 15 1000 400])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');

subplot(1,2,1)
plot(t(1:length(xd)),xd(1,1:length(xd)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(xd)),x(1,1:length(xd)),'-','Color',C21,'LineWidth',lw);hold on;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$x_{1d}$','$x_1$'},'interpreter','latex','fontsize',fontsizeLegend,'Location','best')
subplot(1,2,2)
plot(t(1:length(xd)),xd(2,1:length(xd)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(xd)),x(2,1:length(xd)),'-','Color',C21,'LineWidth',lw);hold on;
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[rad]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$x_{2d}$','$x_2$'},'interpreter','latex','fontsize',fontsizeLegend,'Location','best')

figure('Position', [15 15 1000 400])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');

subplot(1,2,1)
plot(t(1:length(X)),xd(1,1:length(X)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(X)),X(1,1:length(X)),'-','Color',C21,'LineWidth',lw);hold on;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$x_{1d}$','$x_1$'},'interpreter','latex','fontsize',fontsizeLegend,'Location','best')
subplot(1,2,2)
plot(t(1:length(X)),xd(2,1:length(X)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(X)),X(2,1:length(X)),'-','Color',C21,'LineWidth',lw);hold on;
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[rad]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$x_{2d}$','$x_2$'},'interpreter','latex','fontsize',fontsizeLegend,'Location','best')

