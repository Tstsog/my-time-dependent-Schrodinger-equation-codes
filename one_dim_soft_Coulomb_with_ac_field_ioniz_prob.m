% This matlab code solves the time-dependent Schrodinger equation (TDSE) for
% one-dimensional soft-Coulomb potential with an external (AC) field, in which the backward-time centred space (BTCS) difference scheme is used,
% and the ionization probability as function of time is obtained. 
%
% The soft Coulomb potential: V(x) = -1/(1+x^2).
% The time-dependent hamiltonian: H = -0.5*nabla^{2} + 0.5*x^2 + F*x*cos(w*t). 
%
%  Matrix equation in the BTCS scheme:
% (I-lambda*A)*w^(n+1) = w^(n)
%
% An initial condition: psi(x,0) = psi_{0}(x) is ground state wave function of the 1D HO, which is obtained by a finite difference scheme, initially. 
% A boundary condition: psi(-infinity,t) = psi(infinity,t) = 0.
%
% The ionization probability: P_ion = 1 - sum_{bound states}|<\psi_{n}(x)|psi(x,t)>|^{2}.
%
% The atomic unit (au) is used in the calculation. 
% Note that to get more accurate results, one would increase numerical parameters (N, dt, length, etc., )
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% Nov 18, 2023 & University of North Dakota 
%
function [] = one_dim_soft_Coulomb_with_ac_field_ioniz_prob
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
%---
%%% time parameter
%ti = 0.; % t(0)
dt = 0.0250;     % time step
time_iter = 40000;
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
u_mat;
%
V_soft_Coulomb = -1./sqrt(1.+x.^2);
%
H0_ham = -0.5.*u_mat + diag(V_soft_Coulomb);
%
[Vec,En] = eig(H0_ham(2:N,2:N));                                     % Eigenvalue problem
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
[En(1),En(2),En(3),En(4),En(5), En(6), En(7),En(8),En(9),En(10)]
%-0.6698   -0.2748   -0.1466   -0.0580    0.0398    0.1700    0.3223    0.5054    0.7100    0.9433
%           
Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
V1 = Vec(:,1);                         % The unnormalized eigenfunction for the ground state,
%V1 = [0.,;V1,;0.];
n_c = sum(V1.*V1.*dx);
V1 = 1./sqrt(n_c).*V1;
%
V2 = Vec(:,2);                         % The unnormalized eigenfunction for the ground state,
n_c2 = sum(V2.*V2.*dx);
V2 = 1./sqrt(n_c2).*V2;
%
V3 = Vec(:,3);                         % The unnormalized eigenfunction for the ground state,
n_c3 = sum(V3.*V3.*dx);
V3 = 1./sqrt(n_c3).*V3;
%
V4 = Vec(:,4);                         % The unnormalized eigenfunction for the ground state,
n_c4 = sum(V4.*V4.*dx);
V4 = 1./sqrt(n_c4).*V4;
%
V5 = Vec(:,5);                         % The unnormalized eigenfunction for the ground state,
n_c5 = sum(V5.*V5.*dx);
V5 = 1./sqrt(n_c5).*V5;


%%% time dependent calculation begins
%%%
ci = sqrt(-1.);
%En0 = En(1); 
%
% AC field parameters
F_str = 0.100;
omega = 0.150; 
H_pert = F_str.*x;
%
% evolution matrix
unit_I = eye(N+1);
%
psi_t0 = V1;           % initial condition
%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_1 = fopen('one_dim_soft_Coulomb_with_ac_field_ioniz_prob.txt','w');
%%% time loop starts ---
for i = 1:time_iter
    t(i) = i*dt;
%
    psi_old = psi_t0;
    H_ham = -0.5.*u_mat + diag(V_soft_Coulomb) + diag(H_pert).*cos(omega.*t(i));    
    evolution_mat_r = unit_I -  (dt/ci) * H_ham;
    psi_new = evolution_mat_r(2:N,2:N)\psi_old; % % BTCS scheme
    psi_t0 = psi_new;
    %
    N_popul_decay = abs(sum(conj(psi_t0).*V1).*dx).^2 + abs(sum(conj(psi_t0).*V2).*dx).^2 + ...
                    abs(sum(conj(psi_t0).*V3).*dx).^2 + abs(sum(conj(psi_t0).*V4).*dx).^2 + ...
                    abs(sum(conj(psi_t0).*V5).*dx).^2;    
    %
    P_ion_prob = 1. - N_popul_decay;
    %
    output = [i*dt, P_ion_prob ];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.4f\n', output); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_popul_data = fopen('one_dim_soft_Coulomb_with_ac_field_ioniz_prob.txt', 'r');               % 
read_popul_data = textscan(read_popul_data, '%f %f');
ho_time = read_popul_data{1};
Ionization_probability = read_popul_data{2};
%

figure(1)
plot(ho_time, Ionization_probability, 'b', 'LineWidth', 1.5)
xlabel('\mbox{Time} (au)','Interpreter','latex') % ,'fontsize',16
ylabel('Ioniz. prob.','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
%yscale log
box on


%%%
return
end
