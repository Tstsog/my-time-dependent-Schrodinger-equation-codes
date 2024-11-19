% This matlab code solves the time-dependent Schrodinger equation (TDSE) for
% one-dimensional harmonic oscillator (HO) with an external (AC) field, in which the backward-time centred space (BTCS) difference scheme is used,
% and the harmonic generation spectra is obtained. 
%
% The time-dependent hamiltonian: H = -0.5*nabla^{2} + 0.5*x^2 + F*x*cos(w*t). 
%
%  Matrix equation in the BTCS scheme:
% (I-lambda*A)*w^(n+1) = w^(n) + lambda*b^(n+1)
%
% An initial condition: psi(x,0) = psi_{0}(x) is ground state wave function of the 1D HO, which is obtained by a finite difference scheme, initially. 
% A boundary condition: psi(-infinity,t) = psi(infinity,t) = 0.
%
% Dipole moment: d(t) = <psi(x,t)|d|psi(x,t)>
% harmonid generation spectra |d(omega)|^{2} = int_0^time exp(-i*omega*t) * d(d) * dt
%
% The atomic unit (au) is used in the calculation. 
% Note that to get more accurate results, one would increase numerical parameters (N, dt, length, etc., )
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% Nov 18, 2023 & University of North Dakota 
%
function [] = one_dim_har_ocs_with_ac_field_harmonic_generation
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
dt = 0.0100;     % time step
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
V_pot_ho = 0.5.*x.^2;
%
H0_ham = -0.5.*u_mat + diag(V_pot_ho);
%
[Vec,En] = eig(H0_ham(2:N,2:N));                                     % Eigenvalue problem
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
[En(1),En(2),En(3),En(4),En(5), En(6), En(7),En(8),En(9),En(10)]
% 0.499809192293763   1.499045669595092   2.497517893255047   3.495224983449370   4.492166056969852
% 0.5                 1.5                 2.5                 3.5                 4.5               % exact values from E_{n} = (n+1/2), n = 0, 1, 2, ...
%           
Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
V1 = Vec(:,1);                         % The unnormalized eigenfunction for the ground state,
%V1 = [0.,;V1,;0.];
n_c = sum(V1.*V1.*dx);
V1 = 1./sqrt(n_c).*V1;
%

%%% time dependent calculation begins
%%%
ci = sqrt(-1.);
%En0 = En(1); 
%
% AC field parameters
F_str = 0.10;
omega = 0.15; 
H_pert = F_str.*x;
%
% evolution matrix
unit_I = eye(N+1);
%
psi_t0 = V1;           % initial condition
%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_1 = fopen('one_dim_har_ocs_with_ac_field_dipole_ave.txt','w');
%%% time loop starts ---
for i = 1:time_iter
    t(i) = i*dt;
%
    psi_old = psi_t0;
    H_ham = -0.5.*u_mat + diag(V_pot_ho) + diag(H_pert).*cos(omega.*t(i));    
    evolution_mat_r = unit_I -  (dt/ci) * H_ham;
    psi_new = evolution_mat_r(2:N,2:N)\psi_old; % % BTCS scheme
    psi_t0 = psi_new;
    %
    dipole_ave = sum(conj(psi_t0).*x(2:N).*psi_t0.*dx);    
    %
    output = [i*dt, dipole_ave];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.12f\n', output); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_popul_data = fopen('one_dim_har_ocs_with_ac_field_dipole_ave.txt', 'r');               % 
read_popul_data = textscan(read_popul_data, '%f %f');
ho_time = read_popul_data{1};
dipole_ave = read_popul_data{2};
%

figure(1)
plot(ho_time, dipole_ave , 'b', 'LineWidth', 1.5)
xlabel('\mbox{Time} (au)','Interpreter','latex') % ,'fontsize',16
ylabel('Dipole moment','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
%yscale log
box on

%%%%
%t_n = length(ho_time);
%time = zeros(t_n,1);  vals = zeros(t_n,1);
%time_vas(1:t_n) = read_popul_data{1};
%dip_vals(1:t_n) = read_popul_data{2};

%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_2 = fopen('one_dim_har_ocs_with_ac_field_haromic_generation.txt','w');
%
nn = 100.;
%
for ii = 1:nn
    %
    d(ii) = abs(sum(exp(-ci.*omega.*ii.*ho_time).*dipole_ave.*dt)).^2 ;
    %
    output = [ii, d(ii)];
    %
    fprintf(fileID_save_data_2, '%4.4f \t %8.16f\n', output); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_har_gen_data = fopen('one_dim_har_ocs_with_ac_field_haromic_generation.txt', 'r');               % 
read_har_gen_data = textscan(read_har_gen_data, '%f %f');
harmonic_order = read_har_gen_data{1};
harmonic_generation = read_har_gen_data{2};
%

figure(2)
plot(harmonic_order, harmonic_generation , 'b', 'LineWidth', 1.5)
xlabel('harmonic order','Interpreter','latex') % ,'fontsize',16
ylabel('Harmonic generation spectrum','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
yscale log
box on


%%%
return
end
