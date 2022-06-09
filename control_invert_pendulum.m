%% PROGRAM INVERT PENDULUM OVER CAR
%% CLEAN VARIABLES
clc, clear all, close all;

%% SET TIME VARIABLES
to = 0;
t_final = 30;
ts = 0.1;
t = (to:ts:t_final);

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

N = 40; 
%% CREATION EMPTY MATRIX FOR THE EVOLUTION OF THE SYSTEM
x = zeros(4,length(t)+1-N);
x(1,1) = 5; %% INITIAL POSITION OF THE CAR
x(2,1) = 180*(pi/180); %% INITIAL ANGULAR DISPLACEMENT OF THE PENDULUM
x(3,1) = 0; %% INITIAL VELOCITY CAR
x(4,1) = 0; %% INITIAL ANGULAR VELOCITUY OF THE PENDULUM

%% DESIRED STATES
xd = zeros(4,length(t));
xd(1,:) = 0*sin(0.3*t);
xd(2,:) = 0*(pi/180)*sin(0.3*t);

H = x(:,1);
%% INITIALIZATION CONTROL VECTOR
u = zeros(1,length(t));
%% Definicion del horizonte de prediccion


%% Definicion de los limites de las acciondes de control
bounded = [80; -80];
%% Definicion del vectro de control inicial del sistema
vc = zeros(N,1);
H0 = repmat(H,1,N+1)';

%% OPTIMIZATION SOLVER
[f,solver,args] = mpc_invert_pendulum(bounded, N, L, ts);

for k=1:1:length(t)-N
    tic; 
    %% GENERAL VECTOR OF ERROR SYSTEM
    he(:, k) = xd(:,k)-x(:,k);
    
    %% OPTIMAL CONTROLLER SECTION
    %% STATES OF THE OPTIMAL CONTROLLER
    [H0, control] = NMPC_CONTROL(x(:,k), xd(:,:), k, H0, vc, args, solver ,N);

    %% OBTAIN CONTROL VALUES OF THE VECTOR
    u(k) = control(1,1);
   
    %% GET VALUES OF DRONE
    
    x(:,k+1) = f_d(x(:,k),u(:,k),L,ts);
    noise1 = 0.05;
    
    x(:,k+1) = x(:,k+1) + (-noise1+(noise1+noise1)*rand(4,1));
    %% NEW VALUES OPTIMAL CONTROL
    vc = [control(2:end,:);control(end,:)];
    H0 = [H0(2:end,:);H0(end,:)];
    
    %% SAMPLE TIME
    t_sample(k) = toc;
    toc;
end
figure(1)
for k=1:3:length(t)-N
    drawpend(x(:,k),m1,m0,L1);
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
plot(t(1:length(x)),xd(1,1:length(x)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(x)),x(1,1:length(x)),'-','Color',C21,'LineWidth',lw);hold on;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$x_{1d}$','$x_1$'},'interpreter','latex','fontsize',fontsizeLegend,'Location','best')
subplot(1,2,2)
plot(t(1:length(x)),xd(2,1:length(x)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(x)),x(2,1:length(x)),'-','Color',C21,'LineWidth',lw);hold on;
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
plot(t(1:length(x)),xd(1,1:length(x)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(x)),x(1,1:length(x)),'-','Color',C21,'LineWidth',lw);hold on;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[m]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(a)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$x_{1d}$','$x_1$'},'interpreter','latex','fontsize',fontsizeLegend,'Location','best')
subplot(1,2,2)
plot(t(1:length(x)),xd(2,1:length(x)),'-','Color',C19,'LineWidth',lw);hold on;
plot(t(1:length(x)),x(2,1:length(x)),'-','Color',C21,'LineWidth',lw);hold on;
grid minor;

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$[rad]$','interpreter','latex','fontsize',fontsizeLabel)
title({'(b)'},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend({'$x_{2d}$','$x_2$'},'interpreter','latex','fontsize',fontsizeLegend,'Location','best')

