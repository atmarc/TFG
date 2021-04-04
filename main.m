function r = main()
    pot = 1000000;
    for rec=6
        for i=5:6
            N = i * (3^6);
            data = schrodinger_p(N, rec, pot, false);
            disp(data);
            save_to_file([num2str(rec) '_data2.txt'], data);
        end
    end
end

function r = save_to_file(filename, data)
  fileID = fopen(filename,'a');
  fprintf(fileID, [data '\n']);
  fclose(fileID);
end

function r = schrodinger_p(N, rec_lvl, potencial, pbc)
  disp(['Computing: ' num2str(N) ' ' num2str(rec_lvl)]);
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
  
  % Kronecker sum of discrete Laplacians
  L2 = kron(L, I) + kron(I, L);

  % Sierpinski Carpet
  Vext_mat = sierpinski(N, rec_lvl, false, potencial);
  Vext = Vext_mat(:);

  % Hamiltonian
  Hkin = -0.5 * L2;
  Hext = spdiags(Vext, 0, N^2, N^2);
  H = Hkin + Hext; 
  
  % Finding eigenvalues
  precision = 1;
  max_iter = 100000;
  tic
  [PSI,E,ErrorFlag] = lobpcg(randn(N^2, Neig), H, precision, max_iter);
  toc

  % Store image if it is a big calculus
  if rec_lvl > 5
    filename = [num2str(N) '_' num2str(rec_lvl) 'mat.png'];
    imwrite(reshape(PSI, [N,N]), filename);
  end
  
  enregy0 = num2str(E);
  r = [num2str(N) ' ' num2str(rec_lvl) ' ' enregy0 ' ' num2str(ErrorFlag)];
  
end
