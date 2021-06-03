% So random numbers do not repeat
rng('shuffle');

% load values
PSI = load('eigen_mats/vectors_81_1.mat').PSI;
disp('Vectors loaded');
E = load('eigen_mats/values_81_1.mat').E;
disp('Values loaded');

% Number of eigenvalues we are taking, as if we take all of them it gives wierd values
Neig = 1000;
PSI = PSI(:, 1:Neig);
E = E(1:Neig);

Rmax = 0.5;
N = 81;
dx = (Rmax*2)/N;  
x = linspace(-Rmax + dx/2, Rmax - dx/2, N);
[Xmat, Ymat] = meshgrid(x, x); 
X = Xmat(:); Y = Ymat(:);      

% --------------------------------------------------------------------------------------------------------------

sigma = 0.1;
x0 = 0.3; y0 = 0.2;
k = 10;% impulso
A = 0.1;
PSI_init = A * exp(0.5 * (1/sigma.^2) * (- (X - x0).^2 - (Y - y0).^2));
% PSI_init = exp(0.5*(- (X - x0).^2/sigma.^2 - (Y - y0).^2/sigma.^2)) .* exp(sqrt(-1) * k * X);

disp('Vector c computed')

time_values = [0:50] * 0.001;

n_agents = 20;

displacementX = zeros(n_agents, length(time_values));
displacementY = zeros(n_agents, length(time_values));
c = zeros(length(E), 1);

for j=(1:n_agents)
    disp(num2str(j));
    x0 = rand() - Rmax;    
    y0 = rand() - Rmax;    
    PSI_init = A * exp(0.5 * (1/sigma.^2) * (- (X - x0).^2 - (Y - y0).^2));

    for i=1:length(E)  
        c(i) = PSI_init' * PSI(:, i) * dx * dx;
    end

    index = 1;
    for t=time_values
        PSI_t = zeros(N^2, 1);
        for i=1:length(E)  
            PSI_t = PSI_t + c(i) * exp(sqrt(-1) * E(i) * t) * PSI(:, i);
        end
        
        PSI_tn = PSI_t / sqrt(sum(abs(PSI_t).^2 * dx^2));
        
        %xt = sum(X' * PSI_tn.^2 * dx^2);
        
        displacementX(j, index) = sum((X - x0).^2 .* abs(PSI_tn).^2 * dx^2);
        displacementY(j, index) = sum((Y - y0).^2 .* abs(PSI_tn).^2 * dx^2);
        index = index + 1;
        
        % disp(['Time: ' num2str(t)]);
        
        % imagesc(reshape(abs(PSI_t).^2, [N, N]));    
        % pause;

    end
end

disp('Plotting dispacement...');

dispX = zeros(length(time_values), 1);
dispY = zeros(length(time_values), 1);

for j=(1:length(time_values))
    dispX(j) = sum(displacementX(:, j)) / n_agents;
    dispY(j) = sum(displacementY(:, j)) / n_agents;
end

adjust_t = time_values(1:5);
plot(time_values, dispX, time_values, dispY, adjust_t, 0.5 *(sigma^2) + 0.5*(adjust_t.^2)/(sigma^2));

