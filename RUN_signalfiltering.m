% RUN_signalfiltering.m
%
% Computes the resonant frequencies and associated eigenmodes for an array
% of N bubbles graded in size with size factor s using the multipole
% method, and then filters a pure tone signal into a decomposition over the
% eigenmodes.
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

% Define parameters
N = 6;                       % number of domain components / cell bundles
s = 1.05;
for i = 1:N
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
Vol = pi*R.^2;

% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
N_multi = 3;


% Grid for field
gridN = 200;                      % should be increased if N increases
gridMinX1 = floor(-5*R(1));
gridMaxX1 = ceil(cx(end)+6*R(end));
gridMinX2 = -1.1*R(end);
gridMaxX2 = -gridMinX2;
g1 = linspace(gridMinX1, gridMaxX1, gridN);
g2 = [0];
[ g1, g2 ] = meshgrid(g1, g2);
gridPoints = [g1(:) g2(:)]';
gridPointsN = length(gridPoints);
integ_element = (gridMaxX1-gridMinX1)/gridN;


%% Compute initial values
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) min(eig((MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy))));

x = [linspace(0.0001, 0.005, 200), linspace(0.005, 0.025, 300)];
init = [];
for correction = [0 0.0000001i 0.000001i]
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

% clf
% scatter(real(resonances), imag(resonances), 'x')


%% Finding eigenmodes along x2=0
u_store = [];
u_b_store = [];
val_store = [];
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


u = zeros(gridPointsN, 1);
u_b = zeros(gridPointsN, 1);

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
        u_b(j) = I;     % stores the values that are in the bubbles
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

u_store = [u_store, u];
end

J = false(0,0);
numbers = zeros(1,N);
for i = 1:N
    I = (u_b == i);
    numbers(i) = sum(I);
    J = [J, I];
end

gamma = zeros(N);
for n = 1:N
    for m = 1:N
        for k = 1:N
        gamma(n,m) = gamma(n,m) + Vol(k)*dot(u_store(J(:,k),n),u_store(J(:,k),m))/integ_element/numbers(k);
        end
    end
end


norm_u = diag(gamma).^(-0.5);
u_store = u_store*diag(norm_u);
gamma = diag(norm_u)*gamma*diag(norm_u);


%% Response to incoming signal


values_omeg = linspace(0.0001,1.1*real(resonances(end)),300);
norms = [];
alpha_store = [];

for omeg = values_omeg
k = omeg; kb = omeg;
uIn = @(x1, x2) exp(1i*omeg*(x1));  % spatial part of wave only

F = [];
for j = 1:N
    Rad = R(j);
    cenx = cx(j);
    ceny = cy(j);
    for m = -N_multi:N_multi
        fun = @(t) uIn(cenx+Rad*cos(t), ceny+Rad*sin(t)).*exp(-1i*m*t);
        F = [F; integral(fun,-pi,pi)/2/pi];
    end
end

for j = 1:N
    Rad = R(j);
    cenx = cx(j);
    ceny = cy(j);
    for m = -N_multi:N_multi
        fun = @(t) 1i*omeg*cos(t).*uIn(cenx+Rad*cos(t), ceny+Rad*sin(t)).*exp(-1i*m*t);
        F = [F; -delta*integral(fun,-pi,pi)/2/pi];
    end
end
F = real(F);

A = MakeA(R,omeg,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);
V = A\F;

phi = [];
psi = [];
N_terms = 2*N_multi+1;

for i = 1:N
    phi = [phi, V((i-1)*14+1:i*14-7,1)];
    psi = [psi, V((i-1)*14+8:i*14,1)];
end

% Calculate field
green = @(k,x,y) -1i/4*besselh(0,1,k*sqrt(x.^2+y.^2));

u = zeros(gridPointsN, 1);
u_b = zeros(gridPointsN, 4);

parfor j = 1:gridPointsN
    gridPoint = gridPoints(:, j);
    
    % Determine whether we are inside or outside the domains
    I = (gridPoint(1)*ones(1,length(cx))-cx).^2 + (gridPoint(2)*ones(1,length(cy))-cy).^2  <= R.^2 ;
    if sum( I ) > 0
        S = 0;
        I = find(I);
        fun1 = @(t) green(kb, gridPoint(1)-cx(I)-R(I).*cos(t), gridPoint(2)-cy(I)-R(I).*sin(t)).*...
            (exp(-3i*t)*phi(1,I)+exp(-2i*t)*phi(2,I)+exp(-1i*t)*phi(3,I)+phi(4,I)+exp(1i*t)*phi(5,I)+exp(2i*t)*phi(6,I)+exp(3i*t)*phi(7,I));
        S = integral(fun1, -pi, pi);
        u(j) = S;
        u_b(j) = I;     % stores the values that are in the bubbles
    else
        S = 0;
        for i = 1:N
            fun2 = @(t) green(k, gridPoint(1)-cx(i)-R(i).*cos(t), gridPoint(2)-cy(i)-R(i).*sin(t)).*...
                (exp(-3i*t)*psi(1,i)+exp(-2i*t)*psi(2,i)+exp(-1i*t)*psi(3,i)+psi(4,i)+exp(1i*t)*psi(5,i)+exp(2i*t)*psi(6,i)+exp(3i*t)*psi(7,i));
            S = S + integral(fun2, -pi, pi);
        end
        IN = uIn(gridPoint(1),gridPoint(2));
        u(j) = S + IN;
        index(j) = 1;
    end
end


norm_u = 0;
for k = 1:N
    norm_u = norm_u + Vol(k)*dot(u(J(:,k)),u(J(:,k)))/integ_element/numbers(k);
end
norms = [norms norm_u];

RHS = zeros(N,1);
for m = 1:N
    for k = 1:N
        RHS(m) = RHS(m) + Vol(k)*dot(u(J(:,k)),u_store(J(:,k),m))/integ_element/numbers(k);
    end
end

consts = conj(gamma)\RHS;
if size(resonances,1) == 1
    resonances = resonances.';
end
alpha = 2*consts.*imag(resonances).*exp(1i*(real(resonances)-omeg)/2).*sinc((real(resonances)-omeg)/pi);
alpha_store = [alpha_store, alpha];

end



%% Filtering plot

subplot(N+1,1,1)
plot(values_omeg,0.5*log10(norms),'color',[0 0.5 0],'linewidth',1)
ylabel('$\log_{10}\|u\|_X$','interpreter','latex')
set(get(gca,'ylabel'),'rotation',0,'horizontalalignment','right')
xlim([0 values_omeg(end)])
box off
set(gca,'TickLabelInterpreter','latex','xticklabel',[],'FontSize',11)

for n = 1:N
    subplot(N+1,1,n+1)
    plot(values_omeg,abs(alpha_store(n,:)),'linewidth',1)
    xlim([0 values_omeg(end)])
    box off
    set(gca,'TickLabelInterpreter','latex','FontSize',11)
    if n < N
        set(gca,'xticklabel',[])
    else
        xlabel('$\omega_{in}$','interpreter','latex','FontSize',15)
    end    
        
%     ylabel('$|\alpha_n(\Re\omega_n)|$','interpreter','latex')
    ylabel(['$|\alpha_', num2str(n) ,'|$'],'interpreter','latex')

    
    set(get(gca,'ylabel'),'rotation',0,...
        'horizontalalignment','right',...
        'FontSize',15)
end


