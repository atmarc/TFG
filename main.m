function r = main()
  pot = 1000000;
  for N=4:7
    for rec=1:9
      data = schrodinger_p(N*(3^rec), rec, pot, false);
      display([data ' ']);
      %save_to_file('data_zbc_pot_1000000.txt', data);
    end
  end
end

function r = save_to_file(filename, data)
  fileID = fopen(filename,'a');
  fprintf(fileID, [data '\n']);
  fclose(fileID);
end

function r = schrodinger_p(N, rec_lvl, potencial, pbc)
  Neig = 1; % number of eigenvalues to be found
  Rmax = 0.5;

  dx = (Rmax*2)/N;  
  x = linspace(-Rmax + dx/2, Rmax - dx/2, N);   % one-dimensional space lattice
  [Xmat, Ymat] = meshgrid(x, x);                % two-dimensional space lattice
  X = Xmat(:); Y = Ymat(:);                     % all elements of array as a single column
  
  e = ones(N,1);               
  L = spdiags([e -2*e e], -1:1, N, N) / dx^2; % 1D finite difference Laplacian
  I = speye(N);

  if pbc
    % Periodic boundary conditions 
    L(N,1) = 1;
    L(1, N) = 1;
  end
  
  % https://en.wikipedia.org/wiki/Kronecker_sum_of_discrete_Laplacians
  L2 = kron(L, I) + kron(I, L);

  % Sierpinski Carpet
  Vext_mat = sierpinski(N, rec_lvl, false, potencial);
  Vext = Vext_mat(:);

  Hkin = -0.5 * L2;
  Hext = spdiags(Vext, 0, N^2, N^2);
  H = Hkin + Hext;  % Hamiltonian
  
  precision = 1e-3;
  max_iter = 10000;
  tic
  [PSI2,E2,ErrorFlag] = lobpcg(randn(N^2, 2), H, precision, max_iter);
  toc
  display(num2str(E2(1,1)));
  
  
  sigma = 'sa';
  opt.p = 100;
  tic
  [PSI,E] = eigs(H, Neig, sigma, opt);      % Smallest eigenvalue of H
  toc
  display(num2str(E(1,1)));
  
  enregy0 = num2str(E(1,1), 5);
  r = 0;
  %enregy1 = num2str(E(2,1), 5);
  %r = [num2str(N) ' ' num2str(potencial) ' ' num2str(rec_lvl) ' ' enregy0 ' ' enregy1 ' ' num2str(ErrorFlag)];
end