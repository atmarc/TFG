
function r = main()
  pot = 10000;
  for rec=2
      for N=[3^3]
          tic
          data = schrodinger_p(N, rec, 0.5, false, pot);
          disp(data);
          toc
          save_to_file([num2str(rec) '_pbc.txt'], data);
      end
  end
end


function r = schrodinger_p(N, rec_lvl, Rmax, PBC, pot)
disp(['Computing: ' num2str(N) ' ' num2str(rec_lvl)]);
Neig = (N^2 + 1)/2; % number of eigenvalues to be found

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

noise = ((rand(N^2, 1) - 0.5) * 0);
Vext = Vext_mat(:) + noise;
% -------------------------------------

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian

disp('Finding eigenvalues...');
precision = 1e-4;
try
  tic
    [PSI,E,ErrorFlag] = lobpcg(rand(N^2, Neig), H, precision, 10000);
  toc
catch
  tic
  [PSI,E] = eigs(H, Neig, 'sa');
  E = diag(E);
  ErrorFlag = 0;
  toc
end

disp(['Error flag: ' num2str(ErrorFlag)]); % if it doesn't converge with 

disp('Saving matrices...');

save(['data/eigen_mats/vectors_' num2str(N) '_' num2str(rec_lvl) '.mat'], 'PSI');
save(['data/eigen_mats/values_' num2str(N) '_' num2str(rec_lvl) '.mat'], 'E');


enregy0 = num2str(E);
r = [num2str(N) ' ' num2str(rec_lvl) ' ' enregy0 ' ' num2str(ErrorFlag)];

end
