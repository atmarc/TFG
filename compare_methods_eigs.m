%Schrodinger2D.m
disp('Starting program');

% --- Parameters ------------------------------
N = 80;
Neig = 300; % number of eigenvalues to be found
recursion_level = 0;
noise_var = 0;

Rmax = 1/2;
%Rmax = recursion_level / 2;   % --- keep Lmin to 1;
%Rmax = (3^(recursion_level)) / 2;

PBC = true;
% 0 --> no potential | 1 --> Sierpinski carpet | 2 --> chess grid
potential_shape = 0; 
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

%  ---------------- Potential shape ----------------
% Sierpinski Carpet 
if potential_shape == 1
  if PBC
    Vext_mat = sierpinski(N + 1, recursion_level, true);
    Vext_mat = Vext_mat(1:N, 1:N);
  else
    Vext_mat = sierpinski(N, recursion_level, true);
  end
  
  noise = ((rand(N^2, 1) - 0.5) * noise_var);
  Vext = Vext_mat(:) + noise;

% Periodic potential (chess grid)
elseif potential_shape == 2
  if PBC
      Vext_mat = chess_grid(N + 1, recursion_level);
      Vext_mat = Vext_mat(1:N, 1:N);
  else
      Vext_mat = chess_grid(N, recursion_level);
  end
  Vext = Vext_mat(:);

% No potential
else
  Vext_mat = zeros(N,N);
  Vext = Vext_mat(:);
end

% -------------------------------------

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian

disp(['valor N: ' num2str(N)]);
disp('Find eigenvalues...');
% 
% tic
% precision = 1e-3;
% [PSI,E,ErrorFlag] = lobpcg(rand(N^2, Neig), H, precision, 10000);
% toc
% 
tic
%   opt.p = 500; 
%   [PSI,E] = eigs(H, Neig, 'sa', opt);
  [PSI,E] = eigs(H, Neig, 'sa');
  E = diag(E);
toc

% plot([1:Neig], E, 'x');

% E_diff = zeros(length(E) - 1, 1);
% for i=(2:length(E))
%     E_diff(i) = E(i) - E(i - 1);
% end

% plot([1:Neig], E_diff, 'x');

% x = [1:Neig]; plot(x/(2^2), E_2, 'x', x/(4^2), E_4, 'x', x/(8^2), E_8, 'x');
% plot(1:50, E_1, 'x', (1:450)/(3^2), E_2, 'x', (1:4050)/(9^2), E_3, 'x', x/(27^2), E_4, 'x');

%plot(1:50, E_1, x/(3^2), E_2, x/(9^2), E_3, x/(27^2), E_4);

% plot((1:25), E_0, 'o', (1:50)/(2), E_2, 'x', (1:200)/(8), E_4, 'x', (1:800)/(32), E_8, 'x');
% % plot(x/(2^2), E_2, 'x', x/(4^2), E_4, 'x', x/(8^2), E_8, 'x');
% lgnd = legend('no potential', 'iteration 1', 'iteration 2', 'iteration 3'); 
% lgnd.Location = 'southeast'; 
% ylabel('Energy $\left ( \frac{\hbar ^2}{2m} \right )$','Interpreter','latex'); 
% xlabel('Eigenstate','Interpreter','latex'); 
% title('Eigenenergies of the system with periodic potential (pbc)','Interpreter','latex'); 

% xlim([0, Neig]);
% 