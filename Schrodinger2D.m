%Schrodinger2D.m
disp('Starting program');
N = 3^4;
Neig = 5; % number of eigenvalues to be found
Rmax = 0.5;
recursion_level = 4;
PBC = true;

if PBC
  N = N - 1;
end

dx = (Rmax*2)/N;  
x = linspace(-Rmax + dx/2, Rmax - dx/2, N);   % one-dimensional space lattice
[Xmat, Ymat] = meshgrid(x, x);  % two-dimensional space lattice
                                % lattice spacing
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

% https://en.wikipedia.org/wiki/Kronecker_sum_of_discrete_Laplacians
L2 = kron(L, I) + kron(I, L);

% ------ Potencial harmonico ------
%Vext = 0.5*X.^2 + 0.5*Y.^2;

% --------- Sierpinski Carpet ---------
if PBC
  Vext_mat = sierpinski(N + 1, recursion_level, true);
  Vext_mat = Vext_mat(1:N, 1:N);
else
  Vext_mat = sierpinski(N, recursion_level, true);
end

noise = ((rand(N^2, 1) - 0.5) * 0);
Vext = Vext_mat(:) + noise;
%imagesc(reshape(Vext, [N, N]));
%pause;
% -------------------------------------

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian

display('Finding eigenvalues...');

precision = 1e-2;
tic
[PSI,E,ErrorFlag] = lobpcg(rand(N^2, Neig), H, precision, 10000);
toc
display(['Error flag: ' num2str(ErrorFlag)]); % if it doesn't converge with 

max_IPR_wavefun = [];
max_IPR = -1;
min_IPR_wavefun = [];
min_IPR = N;

for i=1:length(diag(E))  
  disp(['Eigenstate ' num2str(i-1) ' energy ' num2str(E(i), 5) '\hbar\omega']); %display result
  PSI_i = PSI(:, i);
  PSI_2 = reshape(PSI_i, [N,N]); 
  PSI_2 = PSI_2 / sign(sum(sum(PSI_2)));
  PSI_2 = PSI_2 / dx;
  
  IPR = sum(sum(sum(PSI_2.^4))*dx^2);
  disp(['IPR: ' num2str(IPR)]);
  data_to_save = [num2str(i-1) ' ' num2str(E(i), 5) ' ' num2str(IPR)];
  
  % save_to_file('data/IPR_data/IPR_data_rec1', data_to_save);
  
  h = pcolor(x,x,PSI_2.^2);
  daspect([1 1 1]);
  colorbar; 
  set(h, 'EdgeColor', 'none');
  pause;
  
  if IPR > max_IPR
    max_IPR = IPR;
    max_IPR_wavefun = PSI_2;
  end
  
  
  if IPR < min_IPR
    min_IPR = IPR;
    min_IPR_wavefun = PSI_2;
  end
  
end

% ------------------------------------------------------
% Display State with max IPR 
disp(['State with max IPR: ' num2str(max_IPR)]);
%imagesc(max_IPR_wavefun.^2); 
%colorbar; 
%pause;

% Display State with min IPR 
disp(['State with min IPR: ' num2str(min_IPR)]);
%imagesc(min_IPR_wavefun.^2); 
%colorbar; 
%pause;
% ------------------------------------------------------

PSI_0 = PSI(:, 1);
E_0 = E(1);
%PSI_0 = PSI_0 * sign(sum(PSI_0));
PSI_mat = reshape(PSI_0, [N,N]); 

for i=floor(N/2 + 1):N
    for j=1:ceil(N/2)
        PSI_mat(j, i) = 0;
    end
end

for i=1:N
    for j=floor(N/2 + 1):N
        PSI_mat(j, i) = 0;
    end
end

%imagesc(PSI_mat);
PSI_c = PSI_mat(:);
%PSI_c = PSI(:, 1) + PSI(:, 2) + PSI(:, 3) + PSI(:, 4);
E_c = (PSI_c' * H * PSI_c) / (PSI_c' * PSI_c);

PSI_prime = H * PSI_c;
E_c_2 = sqrt((PSI_prime' * PSI_prime)/(PSI_c' * PSI_c));
disp(['Energy prime: ' num2str(E_c_2)]);
Eq_18 = (PSI_prime' * PSI_c)^2 / ((PSI_prime' * PSI_prime) * (PSI_c' * PSI_c));
disp(['Comprovar eq 18: ' num2str(Eq_18)]);

disp(['Energy: ' num2str(E_0)]);
disp(['Energy corner: ' num2str(E_c)]);

residual_vec = (PSI_c / sqrt(PSI_c' * PSI_c)) - (PSI_0/norm(PSI_0));
disp(['Residual value: ' num2str(norm(residual_vec))]);







%List of Matlab functions used in the code.

%y = linspace(x1,x2,N) returns N linearly spaced points.
%[X,Y,Z] = meshgrid(xgv,ygv,zgv) produces three-dimensional coordinate arrays. The output coordinate arrays X, Y, and Z contain copies of the grid vectors xgv, ygv, and zgv respectively. The sizes of the output arrays are determined by the length of the grid vectors. For grid vectors xgv, ygv, and zgv of length M, N, and P respectively, X, Y, and Z will have N rows, M columns, and P pages.
%A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
%S = speye(m,n) and S = speye([m n]) form an m-by-n sparse matrix with 1s on the main diagonal.
%K = kron(A,B) returns the Kronecker tensor product of matrices A and B. If A is an m-by-n matrix and B is a p-by-q matrix, then kron(A,B) is an mp-by-nq matrix formed by taking all possible products between the elements of A and the matrix B.
%eigs(A,k,sigma) and eigs(A,B,k,sigma) return k eigenvalues based on sigma, which can take any of the following values… where ‘sa’ means smallest
%disp(X) displays the contents of X without printing the variable name. disp does not display empty variables.
