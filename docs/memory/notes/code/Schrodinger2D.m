%Schrodinger2D.m
Neig = 3; % number of eigenvalues to be found
N = 101;
Rmax = 5.;

x = linspace(-Rmax, Rmax, N);   % one-dimensional space lattice
[Xmat, Ymat] = meshgrid(x, x);        % two-dimensional space lattice
h = x(2) - x(1);                % lattice spacing
dx = h;
X = Xmat(:); Y = Ymat(:);             % all elements of array as a single column

e = ones(N,1);               
L = spdiags([e -2*e e], -1:1, N, N) / h^2; % 1D finite difference Laplacian
I = speye(N);

% https://en.wikipedia.org/wiki/Kronecker_sum_of_discrete_Laplacians
L2 = kron(L, I) + kron(I, L);

Vext = 0.5*X.^2 + 0.5*Y.^2;

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian
[PSI,E] = eigs(H, Neig, 'sa');                      % Smallest eigenvalue of H

for i=1:length(diag(E))
  disp(['Eigenstate ' num2str(i-1) ' energy ' num2str(E(i,i), 5) '\hbar\omega']); %display result
  PSI_2 = reshape(PSI(:, i), [N,N]); h = pcolor(x,x,PSI_2); daspect([1 1 1]); colorbar; set(h, 'EdgeColor', 'none');
  pause 
end

%List of Matlab functions used in the code.

%y = linspace(x1,x2,N) returns N linearly spaced points.
%[X,Y,Z] = meshgrid(xgv,ygv,zgv) produces three-dimensional coordinate arrays. The output coordinate arrays X, Y, and Z contain copies of the grid vectors xgv, ygv, and zgv respectively. The sizes of the output arrays are determined by the length of the grid vectors. For grid vectors xgv, ygv, and zgv of length M, N, and P respectively, X, Y, and Z will have N rows, M columns, and P pages.
%A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
%S = speye(m,n) and S = speye([m n]) form an m-by-n sparse matrix with 1s on the main diagonal.
%K = kron(A,B) returns the Kronecker tensor product of matrices A and B. If A is an m-by-n matrix and B is a p-by-q matrix, then kron(A,B) is an mp-by-nq matrix formed by taking all possible products between the elements of A and the matrix B.
%eigs(A,k,sigma) and eigs(A,B,k,sigma) return k eigenvalues based on sigma, which can take any of the following values… where ‘sa’ means smallest
%disp(X) displays the contents of X without printing the variable name. disp does not display empty variables.