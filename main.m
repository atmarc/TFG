function r = main()
  for N=1:27
    for pot=[1000000]
      for rec=[2]
        data = schrodinger_p(N*(3^rec), rec, pot, false);
        display(data);
        fileID = fopen('data_zbc_pot_1000000.txt','a');
        fprintf(fileID, [data '\n']);
        fclose(fileID);
      end
    end
  end
end

function r = schrodinger_p(N, rec_lvl, potencial, pbc)
  Neig = 2; % number of eigenvalues to be found
  Rmax = 0.5;

  x = linspace(-Rmax, Rmax, N);   % one-dimensional space lattice
  [Xmat, Ymat] = meshgrid(x, x);  % two-dimensional space lattice
  h = x(2) - x(1);                % lattice spacing
  dx = h;
  X = Xmat(:); Y = Ymat(:);       % all elements of array as a single column

  e = ones(N,1);               
  L = spdiags([e -2*e e], -1:1, N, N) / h^2; % 1D finite difference Laplacian
  I = speye(N);

  if pbc
    % Periodic boundary conditions 
    L(N,1) = 1;
    L(1, N) = 1;
  end
  
  % https://en.wikipedia.org/wiki/Kronecker_sum_of_discrete_Laplacians
  L2 = kron(L, I) + kron(I, L);

  % --------- Sierpinski Carpet ---------
  %%{
  Vext_mat = sierpinski(N, rec_lvl, false, potencial);
  Vext = Vext_mat(:);
  %%}
  % -------------------------------------

  Hkin = -0.5 * L2;
  Hext = spdiags(Vext, 0, N^2, N^2);
  H = Hkin + Hext;  % Hamiltonian

  
  sigma = 'sa';
  opt.p = 200;
  %display('Looking for eigenvalues...');
  [PSI,E] = eigs(H, Neig, sigma, opt);      % Smallest eigenvalue of H

  enregy0 = num2str(E(1,1), 5);
  enregy1 = num2str(E(2,2), 5);
  
  r = [num2str(N) ' ' num2str(potencial) ' ' num2str(rec_lvl) ' ' enregy0 ' ' enregy1];
end