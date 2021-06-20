% So random numbers do not repeat
rng('shuffle');

% load values
iteration = 0;
folder = 'zbc';
PSI = load(['data/eigen_mats/' folder '/vectors_81_' num2str(iteration) '.mat']).PSI;
disp('Vectors loaded');
E = load(['data/eigen_mats/' folder '/values_81_' num2str(iteration) '.mat']).E;
disp('Values loaded');

% Number of eigenvalues we are taking, as if we take all of them it gives wierd values
Neig = 1000;
PSI = PSI(:, 1:Neig);
E = E(1:Neig);
Rmax = (3^(iteration)) / 2;
N = 81;
dx = (Rmax*2)/N;  
x = linspace(-Rmax + dx/2, Rmax - dx/2, N);
[Xmat, Ymat] = meshgrid(x, x); 
X = Xmat(:); Y = Ymat(:);      

% --------------------------------------------------------------------------------------------------------------

sigma = 0.1;
x0 = 0; y0 = 0;
k = 10;% impulso
A = 1;
PSI_init = A * exp(0.5 * (1/sigma.^2) * (- (X - x0).^2 - (Y - y0).^2));
% PSI_init = exp(0.5*(- (X - x0).^2/sigma.^2 - (Y - y0).^2/sigma.^2)) .* exp(sqrt(-1) * k * X);

time_values = (0:100) * 0.0002;

n_agents = 1;
disp(   ['Computing ' num2str(n_agents) ' agents']);

displacementX = zeros(n_agents, length(time_values));
displacementY = zeros(n_agents, length(time_values));
c = zeros(length(E), 1);

for j=(1:n_agents)
    disp(num2str(j));
%      x0 = rand()*Rmax*2 - Rmax;    
%      y0 = rand()*Rmax*2 - Rmax;    
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
        
        displacementX(j, index) = sum((X - x0).^2 .* abs(PSI_tn).^2 * dx^2);
        displacementY(j, index) = sum((Y - y0).^2 .* abs(PSI_tn).^2 * dx^2);
        index = index + 1;
        
        disp(['Time: ' num2str(t)]);
        imagesc(reshape(abs(PSI_t).^2, [N, N]));    
        pause;
        
%         h = pcolor(x,x,PSI_2);
%       daspect([1 1 1]);
%       colorbar; 
%       set(h, 'EdgeColor', 'none');
%       pause;

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

% for j=(1:length(time_values))
%     save_to_file('data/time_evolution_disp/disp_0_smallt.txt', [num2str(dispX(j)) ' ' num2str(dispY(j))]);
% end

