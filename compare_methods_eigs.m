%Schrodinger2D.m
disp('Starting program');

% --- Parameters ------------------------------
N = 3^4;
Neig = N; % number of eigenvalues to be found
recursion_level = 3;
Rmax = 1 / 2;
noise_var = 0;
%Rmax = (3^(recursion_level)) / 2   % --- keep Lmin to 1;
PBC = false;
% ---------------------------------------------

if PBC
  N = N - 1;
end

dx = (Rmax*2)/N;  
x = linspace(-Rmax + dx/2, Rmax - dx/2, N);   % one-dimensional space lattice
[Xmat, Ymat] = meshgrid(x, x);  % two-dimensional space lattice
h = dx;
X = Xmat(:); Y = Ymat(:);       % all elements of array as a single column

e = ones(N,1);

L = spdiags([e -2*e e], -1:1, N, N);

% Periodic boundary conditions 
if PBC
  L(N,1) = 1;
  L(1, N) = 1;
end

L = L / h^2; % 1D finite difference Laplacian
I = speye(N);
L2 = kron(L, I) + kron(I, L);

% --------- Sierpinski Carpet ---------
if PBC
  Vext_mat = sierpinski(N + 1, recursion_level, true);
  Vext_mat = Vext_mat(1:N, 1:N);
else
  Vext_mat = sierpinski(N, recursion_level, true);
end

noise = ((rand(N^2, 1) - 0.5) * noise_var);
Vext = Vext_mat(:) + noise;
%imagesc(reshape(Vext, [N, N]));
%pause;
% -------------------------------------

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian

disp(['valor N: ' num2str(N)]);
disp('Find eigenvalues...');

% tic
% precision = 1e-3;
% [PSI,E,ErrorFlag] = lobpcg(rand(N^2, Neig), H, precision, 10000);
% toc

tic
  opt.p = 6561;
  %[PSI,E] = eigs(H, Neig, 'sa', opt);
  [PSI,E] = eigs(H, Neig, 'sa', opt);
  E = diag(E);
toc

plot([1:Neig], E, 'x');