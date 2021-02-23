n = 4;
M = zeros(n,n);
I = speye(n);
for i=1:n
  for j=1:n
    M(i,j) = (j*10) + i;
  endfor
endfor

display(M);
K1 = kron(I, M);
K2 = kron(M, I);
display(full(K1));
display('-----------------------------');
display(full(K2));

F = zeros(n^2,n^2);
for i=1:(n*n)
  for j=1:(n*n)
    F(i,j) = K1(i,j)*100 + K2(i,j);
  endfor
endfor
display('-----------------------------');
display(F);