%Schrodinger2D.m
disp('Starting program');

% --- Parameters ------------------------------
N = 3^3;
Neig = N;%(N^2 + 1)/2; % number of eigenvalues to be found
recursion_level = 1;
Rmax = 1 / 2;
noise_var = 0;
% Rmax = (3^(recursion_level)) / 2   % --- keep Lmin to 1;
PBC = true;
% ---------------------------------------------

if PBC
    N = N - 1;
end

dx = (Rmax*2)/N;  
x = linspace(-Rmax + dx/2, Rmax - dx/2, N);   % one-dimensional space lattice
[Xmat, Ymat] = meshgrid(x, x);  % two-dimensional space lattice
h = dx;
X = Xmat(:); 
Y = Ymat(:);

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
    Vext_mat = sierpinski(N + 1, recursion_level, true);
    Vext_mat = Vext_mat(1:N, 1:N);
else
    Vext_mat = sierpinski(N, recursion_level, true);
end

noise = ((rand(N^2, 1) - 0.5) * noise_var);
Vext = Vext_mat(:) + noise;
%imagesc(reshape(Vext, [N, N]));
%pause;
% -------------------------------------

Hkin = -0.5 * L2;
Hext = spdiags(Vext, 0, N^2, N^2);
H = Hkin + Hext;  % Hamiltonian

disp('Finding eigenvalues...');

precision = 1e-2;
tic
    %[PSI,E,ErrorFlag] = lobpcg(rand(N^2, Neig), H, precision, 10000);
    [PSI,E] = eigs(H, Neig, 'sa');
    E = diag(E);
toc

% --------------------------------------------------------------------------------------------------------------

c = zeros(length(E), 1);

sigma = 0.2;
x0 = 0; y0 = 0;
k = 10;% impulso
% PSI_init = exp(0.5*(- (X - x0).^2/sigma.^2 - (Y - y0).^2/sigma.^2)) .* exp(sqrt(-1) * k * X);
%PSI_init = PSI(:, 1);
sigma = 0.1;
PSI_init = exp(0.5 * (1/sigma.^2) * (- (X - x0).^2 - (Y - y0).^2));

for i=1:length(E)  
    disp(['Eigenstate: ' num2str(i-1) ' Energy: ' num2str(E(i), 5)]); %display result
    PSI_i = PSI(:, i);

    %PSI_i = PSI_i / sqrt(sum(PSI_i.^2)*dx*dx);
    %PSI_i = PSI_i / sign(sum(PSI_i));

    c(i) = PSI_init' * PSI_i * dx * dx;
end



time_values = [0:50] * 0.005;
displacementX = zeros(1, length(time_values));
displacementY = zeros(1, length(time_values));
index = 1;
for t=time_values
    PSI_t = zeros(N^2, 1);
    for i=1:length(E)  
        PSI_t = PSI_t + c(i) * exp(sqrt(-1) * E(i) * t) * PSI(:, i);
    end
    
    PSI_tn = PSI_t / sqrt(sum(abs(PSI_t).^2 * dx^2));
    
    %xt = sum(X' * PSI_tn.^2 * dx^2);
    
    displacementX(1, index) = sum((X - x0).^2 .* abs(PSI_tn).^2 * dx^2);
    displacementY(1, index) = sum((Y - y0).^2 .* abs(PSI_tn).^2 * dx^2);
    index = index + 1;
    
    disp(['Time: ' num2str(t)]);
    
    imagesc(reshape(abs(PSI_t).^2, [N, N]));
    pause;
    
    %drawnow 
    % Capture the plot as an image 
    %frame = getframe(my_figure); 
    %im = frame2im(frame); 
    %[imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    %if t == 0
    %    imwrite(imind,cm,filename,'gif','DelayTime',0.1,'Loopcount',inf); 
    %else 
    %    imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    %end 
    
end

disp('Plotting dispacement...');

adjust_t = time_values(1:5);

plot(time_values, displacementX, time_values, displacementY, adjust_t, 0.5 *(sigma^2) + 0.5*(adjust_t.^2)/(sigma^2));



