% RUN_travellingwave.m
%
% Computes the resonant frequencies and associated eigenmodes for an array
% of N bubbles graded in size with size factor s using the multipole
% method, and then...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B
%
% Used to create Fig 5 of:
% A fully-coupled subwavelength resonance approach to modelling the passive
% cochlea (2019)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

%% Define parameters

N = 50;                       % number of domain components / cell bundles
s = 1.05;
R(1) = 1;
for i = 2:N
    R(i) = s^(i-1);
end


% Material parameters
% rho_b, kappa_b : resonators
% rho0, kappa0 : background
high=7000;         
   
rho_b=1;
rho0=high*1;
kappa_b=1;
kappa0=high*1;

% High contrast parameter \delta
delta=rho_b/rho0;

cx = R(1)*ones(1,N);
if N > 1
    for i = 2:N
        cx(i) = cx(i-1) + 2*R(i-1) + R(i);
    end
end
cy = zeros(1,N);

% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
N_multi = 3;

%% Compute initial conditions
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) min(eig((MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy))));

x = [linspace(0.0001, 0.005, 200), linspace(0.005, 0.025, 400)];
init = [];
for correction = [0 0.0000001i 0.000001i 0.0001i]
    y = zeros(1, length(x));
    for i = 1:length(x)
        y(i) = abs(f(x(i) - 0.00000001i/x(i) - correction));
    end
    for i = 2:length(x)-1
        if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e-8)
            init = [init x(i)];
        end
    end
end

if length(init) < length(R)
    disp('WARNING: fewer than N initial guesses created')
end

init = sort(init);

%% Use Muller's method to find the resonances

distTol = 5e-5; fTol = 1e-5; iterMax = 10;
resonances = [];
n = 1;
for initGuess = init
        
    z0 = initGuess;
    z1 = initGuess - 0.00001i;
    z2 = initGuess - 0.00002i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e-7
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi \n'], real(res), imag(res))
       resonances = [resonances res];
       n = n + 1;
    end
end

% scatter(real(resonances), imag(resonances), 'x')

%% Computing approx. value of eigenmodes on each bubble
u_store = [];
Vol = pi*R.^2;
N_res = length(resonances);

for m = 1:N_res
omega = resonances(m);
k = omega;                          % assumed v = 1
kb = omega;                         % assumed v_b = 1

A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);
[V, D] = eig(A);
[~,permutation] = sort(diag(D));
D=D(permutation,permutation);V = V(:,permutation);

phi = [];
psi = [];
N_terms = 2*N_multi+1;

for i = 1:N
    phi = [phi, V((i-1)*14+1:i*14-7,1)];
    psi = [psi, V((i-1)*14+8:i*14,1)];
end

% Calculate field
green = @(k,x,y) -1i/4*besselh(0,1,k*sqrt(x.^2+y.^2));

gridPoints = [cx; cy];
gridPointsN = length(gridPoints);

u = zeros(gridPointsN, 1);
u_b = zeros(gridPointsN, 4);

parfor j = 1:gridPointsN
    gridPoint = gridPoints(:, j);
    
    % Determine whether we are inside or outside the domains
    I = (gridPoint(1)*ones(1,length(cx))-cx).^2 + (gridPoint(2)*ones(1,length(cy))).^2  <= R.^2 ;
    if sum( I ) > 0
        S = 0;
        I = find(I);
        fun1 = @(t) green(kb, gridPoint(1)-cx(I)-R(I).*cos(t), gridPoint(2)-cy(I)-R(I).*sin(t)).*...
            (exp(-3i*t)*phi(1,I)+exp(-2i*t)*phi(2,I)+exp(-1i*t)*phi(3,I)+phi(4,I)+exp(1i*t)*phi(5,I)+exp(2i*t)*phi(6,I)+exp(3i*t)*phi(7,I));
        S = integral(fun1, -pi, pi);
        u(j) = S;
        u_b(j,:) = [gridPoint(1), gridPoint(2), S, I];     % stores the values that are in the bubbles
    else
        S = 0;
        for i = 1:N
            fun2 = @(t) green(k, gridPoint(1)-cx(i)-R(i).*cos(t), gridPoint(2)-cy(i)-R(i).*sin(t)).*...
                (exp(-3i*t)*psi(1,i)+exp(-2i*t)*psi(2,i)+exp(-1i*t)*psi(3,i)+psi(4,i)+exp(1i*t)*psi(5,i)+exp(2i*t)*psi(6,i)+exp(3i*t)*psi(7,i));
            S = S + integral(fun2, -pi, pi);
        end
        u(j) = S;
        index(j) = 1;
    end
end

% normalisation
avg_val = zeros(N,1);
for i = 1:N
    I = (u_b(:,4) == i);
    total_vals = sum(u_b(I,3));
    number = sum(I);
    avg_val(i) = total_vals/number;
end

norm_u = dot(Vol,avg_val.*conj(avg_val));
norm_u = sqrt(norm_u);

u = u/norm_u;
u_store = [u_store, u];
end



%% Response to initial data
% Note that the approx values of the eigenmodes on each bubble are stored 
% as the columns of u_store

u_0 = 0;
v_0 = 0.05;

for n = 1:N_res
    beta(n) = Vol*u_0*u_store(:,n);
end

for n = 1:N_res
   alpha(n) = Vol*v_0*u_store(:,n)/real(resonances(n)) - ...
       imag(resonances(n))*beta(n)/real(resonances(n)); 
end

u_store1dim = [];
for m = 1:N_res
omega = resonances(m);
k = omega;                          % assumed v = 1
kb = omega;                         % assumed v_b = 1

A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);
[V, D] = eig(A);
[~,permutation] = sort(diag(D));
D=D(permutation,permutation);V = V(:,permutation);

phi = [];
psi = [];
N_terms = 2*N_multi+1;

for i = 1:N
    phi = [phi, V((i-1)*14+1:i*14-7,1)];
    psi = [psi, V((i-1)*14+8:i*14,1)];
end

% Calculate field
green = @(k,x,y) -1i/4*besselh(0,1,k*sqrt(x.^2+y.^2));

% Grid for field
gridN = 500;
gridMinX1 = -R(1);
gridMaxX1 = cx(end) + 2*R(end);
g1 = linspace(gridMinX1, gridMaxX1, gridN);
g2 = zeros(1,gridN); 
gridPoints = [g1; g2];
gridPointsN = length(gridPoints);

u = zeros(gridPointsN, 1);

parfor j = 1:gridPointsN
    gridPoint = gridPoints(:, j);
    
    % Determine whether we are inside or outside the domains
    I = (gridPoint(1)*ones(1,length(cx))-cx).^2 + (gridPoint(2)*ones(1,length(cy))).^2  <= R.^2 ;
    if sum( I ) > 0
        S = 0;
        I = find(I);
        fun1 = @(t) green(kb, gridPoint(1)-cx(I)-R(I).*cos(t), gridPoint(2)-cy(I)-R(I).*sin(t)).*...
            (exp(-3i*t)*phi(1,I)+exp(-2i*t)*phi(2,I)+exp(-1i*t)*phi(3,I)+phi(4,I)+exp(1i*t)*phi(5,I)+exp(2i*t)*phi(6,I)+exp(3i*t)*phi(7,I));
        S = integral(fun1, -pi, pi);
        u(j) = S;
    else
        S = 0;
        for i = 1:N
            fun2 = @(t) green(k, gridPoint(1)-cx(i)-R(i).*cos(t), gridPoint(2)-cy(i)-R(i).*sin(t)).*...
                (exp(-3i*t)*psi(1,i)+exp(-2i*t)*psi(2,i)+exp(-1i*t)*psi(3,i)+psi(4,i)+exp(1i*t)*psi(5,i)+exp(2i*t)*psi(6,i)+exp(3i*t)*psi(7,i));
            S = S + integral(fun2, -pi, pi);
        end
        u(j) = S;
        index(j) = 1;
    end
end

u_store1dim = [u_store1dim, real(u)];

end


%%

times = 1000:2000:9000;
p = 1;
for t = times
    u = 0;
    for n = 2:N_res
        u = u + u_store1dim(:,n)*exp(imag(resonances(n))*t)*(alpha(n)*sin(real(resonances(n))*t)+beta(n)*cos(real(resonances(n))*t));
    end
    subplot(length(times),1,p)
    plot(g1,real(u));
    title(['t = ' num2str(t)])
    ylim([-500 500]), 
    xlim([gridMinX1 gridMaxX1])
    p = p + 1;
end
