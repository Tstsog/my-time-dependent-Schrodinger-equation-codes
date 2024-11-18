% This matlab code solves the time-dependent Schrodinger equation (TDSE) for
% one-dimensional harmonic oscillator (HO), in which the backward-time centred space (BTCS) difference scheme is used.  
%
%  Matrix equation in the BTCS scheme:
% (I-lambda*A)*w^(n+1) = w^(n) + lambda*b^(n+1)
%
% An initial condition: psi(x,0) = psi_{0}(x) is ground state wave function of the 1D HO, which is obtained by a finite difference scheme, initially. 
% A boundary condition: psi(-infinity,t) = psi(infinity,t) = 0.
%
% The time-independent and time-dependent problems have analytical solutions, with which an numerical solutions are compared, as well.   
% The atomic unit (au) is used in the calculation.
% Note that to get more accurate results, one would increase numerical parameters (N, dt, length, etc., )
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% Nov 18, 2023 & University of North Dakota 
%
function [] = one_dim_har_ocs_tdse_n0
clc; 
format short
%
% grid and initial condition
a = -10.; % x(0)
b = 10.; % x(N+1)
N = 256;  % number of grid point of x axis
dx = (b-a)/N; % step size in x
% ---
x = zeros(N+1,1); % total number of points is N+1
for i = 1:N+1
    x(i) = a + (i-1)*dx;
end
x;
% ---
psi_n0_exact = (1/pi).^(1/4).*exp(-0.5.*x.^2);
%---
%%% time parameter
%ti = 0.; % t(0)
dt = 0.01;     % time step
time_iter = 50;
%
lambda = 1/dx^2; 
%  matrix equation is
% (I-lambda*A)*w^(n+1) = w^(n) + lambda*b^(n+1)
u_mat = zeros(N+1,N+1);
for i = 2:N
    u_mat(i,i-1) = lambda;
    u_mat(i,i) = -2.*lambda;
    u_mat(i,i+1) = lambda;
end
u_mat(1,1) = -2.*lambda; u_mat(N+1,N+1) = -2.*lambda ; 
u_mat(1,2) = lambda; u_mat(N+1,N) = lambda;
%u_mat;
%
V_pot_ho = 0.5.*x.^2;
%
H0_ham = -0.5.*u_mat + diag(V_pot_ho);
%
[Vec,En] = eig(H0_ham(2:N,2:N));                                     % Eigenvalue problem
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
[En(1),En(2),En(3),En(4),En(5)]
% 0.499809192293763   1.499045669595092   2.497517893255047   3.495224983449370   4.492166056969852
% 0.5                 1.5                 2.5                 3.5                 4.5               % exact values from E_{n} = (n+1/2), n = 0, 1, 2, ...
%           
Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
V1 = Vec(:,1);                         % The unnormalized eigenfunction for the ground state,
%V1 = [0.,;V1,;0.];
n_c = sum(V1.*V1.*dx);
V1 = 1./sqrt(n_c).*V1;
%

figure(1)
hold on
plot(x(2:N), V1, 'b', 'LineWidth',1.5) % numerical wave function for ground state
plot(x(2:N), psi_n0_exact(2:N), 'ro')  % analytical wave function 
hold off
xlabel('x\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$\psi_{0}(x)$','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on

%%% time dependent calculation begins
% exact soluton for time-dependent Schrodinger equation with 
% psi(x,0) = psi_{0}(x) initial condition is 
% psi(x,t) = (1/pi)^(1/4)*exp(-i*E_{0}*t)*exp(-0.5*x^2);
%
%%%
ci = sqrt(-1.);
En0 = En(1); 
% evolution matrix
unit_I = eye(N+1);
evolution_mat_r = unit_I -  (dt/ci) * H0_ham;
%
psi_t0 = V1;           % initial condition
%%%%%%%%%%%%%%%%%%%%%%%%
%%% time loop starts ---
for i = 1:time_iter
%    t(i) = i*dt;
%
    psi_old = psi_t0;
    psi_new = evolution_mat_r(2:N,2:N)\psi_old; % % BTCS scheme
    psi_t0 = psi_new;
    %
end

%%%
figure(2)
hold on 
%plot(x(2:N), V1, 'b', 'LineWidth',1.5) % numerical wave function for ground state
%plot(x(2:N), psi_n0_exact(2:N), 'go')  % analytical wave function 
%
plot(x(2:N), conj(psi_t0).*psi_t0, 'b', 'LineWidth',1.5)
plot(x(2:N), conj((1./pi).^(1/4).*exp(-ci.*En0.*time_iter).*exp(-0.5.*x(2:N).^2)).*((1./pi).^(1/4).*exp(-ci.*En0.*time_iter).*exp(-0.5.*x(2:N).^2)), 'ro')
xlabel('x\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$|\psi_{0}(x,t)|^{2}$','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on



%%%
return
end
