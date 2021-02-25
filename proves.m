N = 10000;

% Random symetric matrix
d = 1000000*rand(N,1);
t = triu(bsxfun(@min,d,d.').*rand(N),1); 
M = diag(d)+t+t.';

tic
[PSI,E] = eigs(M, 2, 'sa');
toc
tic
[PSI2,E2,ErrorFlag] = lobpcg(rand(N, 2), M, 1e-4, 10000);
toc
display([num2str(E(1,1)) ' ' num2str(E2(1,1))])