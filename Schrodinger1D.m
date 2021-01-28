%Schrodinger1D.m
N = 101; % number of points
Rmax = 6.; % maximal distance

Neig = 3; % number of eigenvalues to be found

x = linspace(-Rmax, Rmax, N);   % one-dimensional space lattice
h = x(2) - x(1);                % lattice spacing

% kinetic energy
e = ones(N,1);               
L = spdiags([e -2*e e], -1:1, N, N) / h^2; % 1D finite difference Laplacian
Hkin = -0.5 * L;

% extrenal potential
Vext = 0.5 * x'.^2;
Hext = spdiags(Vext, 0, N, N);

H = Hkin + Hext;

% Smallest eigenvalues of H
[PSI,E] = eigs(H, Neig, 'sa');

% ground state
psi = PSI(:,1); % wave function
psi = psi / sign(sum(psi)); % change sign if needed
Etot = psi' * H * psi ./ (psi' * psi) 
Ekin = psi' * Hkin * psi ./ (psi' * psi)
Eext  = psi' * Hext * psi ./ (psi' * psi)
plot(x, psi)

% first excited state
psi = PSI(:,2); % wave function
psi = psi / sign(sum(psi)); % change sign if needed
Etot = psi' * H * psi ./ (psi' * psi) 
Ekin = psi' * Hkin * psi ./ (psi' * psi)
Eext  = psi' * Hext * psi ./ (psi' * psi)

hold
plot(x, psi)

