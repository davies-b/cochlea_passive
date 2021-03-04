% RUN_tonotopic.m
%
% Computes the resonant frequencies and associated eigenmodes for an array
% of N bubbles graded in size with size factor s using the multipole
% method, then studies the locations of maximum amplitude of each mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B
%
% Used to create Fig 4 of:
% A fully-coupled subwavelength resonance approach to modelling the passive
% cochlea (2019)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

%% Define parameters

N = 50;                       % number of domain components / cells bundles
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

% High contrast parameters \delta
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
%     pause
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


%% Computing eigenmodes over the plane
u_store = [];
U_store = [];
norm_store = [];
Vol = pi*R.^2;
N_res = length(resonances);

for m = 1:N_res
    m
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
gridN = N*10+1;            
gridMinX1 = floor(-10-R(1));
gridMaxX1 = ceil(cx(end)+R(end)+10);
gridMinX2 = -5;
gridMaxX2 = -gridMinX2;
g1 = linspace(gridMinX1, gridMaxX1, gridN);
% g2 = linspace(gridMinX2, gridMaxX2, gridN);
g2 = zeros(1,gridN);
gridPoints = [g1(:) g2(:)]';

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
norm_store = [norm_store, norm_u];
U_store = [U_store, avg_val];

figure
plot(g1, abs(u)), xlim([gridMinX1, gridMaxX1])
end


%% Tonotopic map investigation

% ignore the lower freq. modes, for which the trend does not hold
start = 18;             % 18 corresponds to N = 50

pos = zeros(1,N_res);
loc = pos;
for n = 1:N_res
    [~,pos(n)] = max(abs(U_store(:,n)));
    loc(n) = cx(pos(n));
end

freqs = real(resonances(start:end));
f = @(b,x) b(1).*exp(b(2).*x)+b(3);                                     
B = fsolve(@(b) norm(freqs - f(b,loc(start:end))), [0.01; -0.005; 0.01])   

figure
plot(loc(start:end), freqs, 'x')
hold on
plot(loc(start:end), f(B,loc(start:end)), '-r')
plot(loc(1:start-1),real(resonances(1:start-1)),'x')
text(200, 0.016, sprintf('%.4f\\cdote^{%.4f\\cdotx}%+.4f', B),'color','red')
