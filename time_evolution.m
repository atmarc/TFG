%Schrodinger2D.m
disp('Starting program');

% --- Parameters ------------------------------
N = 3^4;
Neig = N; % number of eigenvalues to be found
recursion_level = 4;
Rmax = 1 / 2;
noise_var = 0;
% Rmax = (3^(recursion_level)) / 2   % --- keep Lmin to 1;
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

display('Finding eigenvalues...');

precision = 1e-4;
tic
    [PSI,E,ErrorFlag] = lobpcg(rand(N^2, Neig), H, precision, 10000);
toc

display(['Error flag: ' num2str(ErrorFlag)]); % if it doesn't converge with 

c = zeros(length(E), 1);

PSI_init = ...

for i=1:length(E)  
    disp(['Eigenstate: ' num2str(i-1) ' Energy: ' num2str(E(i), 5)]); %display result
    PSI_i = PSI(:, i);

    PSI_i = PSI_i / sqrt(sum(PSI_i.^2)*dx*dx);
    PSI_i = PSI_i / sign(sum(PSI_i));

    c(i) = PSI_init' * PSI_i * dx*dx;
end


for t=1:100
    PSI_t = zeros(length(E), 1);
    for i=1:length(E)  
        PSI_t = PSI_t + c(i) * exp(sqrt(-1) * E(i) * t) * PSI(:, i);
    end

    imagesc(reshape(PSI_t, [N, N]));
    pause;
    
end
