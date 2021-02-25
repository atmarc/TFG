import numpy as np

Neig = 1 # number of eigenvalues to be found
N = 3^5
Rmax = 0.5
recursion_level = 5

x = np.linspace(-Rmax, Rmax, N) # one-dimensional space lattice
Xmat, Ymat = np.meshgrid(x, x)  # two-dimensional space lattice
dx = x[2] - x[1]                # lattice spacing
h = dx
X = Xmat.flatten('F')
Y = Ymat.flatten('F')           # all elements of array as a single column

e = ones(N,1);

L = spdiags([e -2*e e], -1:1, N, N);

% Periodic boundary conditions 
L(N,1) = 1;
L(1, N) = 1;

L = L / h^2; % 1D finite difference Laplacian

I = speye(N);

L2 = kron(L, I) + kron(I, L);

% ------ Potencial harmonico ------
%Vext = 0.5*X.^2 + 0.5*Y.^2;


% --------- Sierpinski Carpet ---------
%%{
Vext_mat = sierpinski(N, recursion_level, false);
imshow(Vext_mat);
pause;
Vext = Vext_mat(:);
%%}
% -------------------------------------

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian

display('Finding eigenvalues...');

opt.p = 100;
%sigma = 'si'; 
sigma = 'sa'; 
[PSI,E] = eigs(H, Neig, sigma, opt);      % Smallest eigenvalue of H

for i=1:length(diag(E))
  disp(['Eigenstate ' num2str(i-1) ' energy ' num2str(E(i,i), 5) '\hbar\omega']); %display result
  comp_value = 1/(1/3^recursion_level)^2;
  %display([num2str(comp_value, 5) ' ' num2str(E(i,i),5)]);
  PSI_2 = reshape(PSI(:, i), [N,N]); 
  PSI_2 = PSI_2 / sign(sum(sum(PSI_2)));
  PSI_2 = PSI_2 / dx;
  h = pcolor(x,x,PSI_2);
  %sum(sum((((Vext_mat > 0).*PSI_2).^2)))*dx*dx
  %sum(sum((((Vext_mat > -1).*PSI_2).^2)))*dx*dx
  
  daspect([1 1 1]);
  colorbar; 
  set(h, 'EdgeColor', 'none');
  pause;
end