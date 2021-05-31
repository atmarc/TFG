
function r = main()
  pot = 10000;
  for rec=6
      for N=[3^6 3^6*2 3^6*3 3^6*4 3^6*5 3^6*6 3^6*7 3^6*8]
          tic
          data = schrodinger_p(N, rec, false, pot);
          disp(data);
          toc
          save_to_file(['rmax_proves_' num2str(rec) '_pbc.txt'], data);
      end
  end
end


function r = schrodinger_p(N, rec_lvl, PBC, pot)
disp(['Computing: ' num2str(N) ' ' num2str(rec_lvl)]);
Neig = 1; % number of eigenvalues to be found
Rmax = (3^(recursion_level)) / 2;

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
  Vext_mat = sierpinski(N + 1, rec_lvl, false, pot);
  Vext_mat = Vext_mat(1:N, 1:N);
else
  Vext_mat = sierpinski(N, rec_lvl, false, pot);
end

Vext = Vext_mat(:);
% -------------------------------------

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian

disp('Finding eigenvalues...');
precision = 1e-3;
try
  [PSI,E,ErrorFlag] = lobpcg(rand(N^2, Neig), H, precision, 10000);
catch
  disp(['Computing values with eigs...']);  
  [PSI,E] = eigs(H, Neig, 'sa');
  E = diag(E);
  ErrorFlag = 0;
end

disp(['Error flag: ' num2str(ErrorFlag)]); % if it doesn't converge with 

disp('Saving matrices...');


enregy0 = num2str(E(1));
r = [num2str(N) ' ' num2str(rec_lvl) ' ' enregy0 ' ' num2str(ErrorFlag)];

end
