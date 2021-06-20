Neig = 300;

% Values computed with Rmax = 3^(recursion_level) / 2;
% E_0 = load('data/eigen_mats/pbc_L_variable/values_81_0.mat').E(1:Neig);
% E_1 = load('data/eigen_mats/pbc_L_variable/values_81_1.mat').E(1:Neig);
% E_2 = load('data/eigen_mats/pbc_L_variable/values_81_2.mat').E(1:Neig);
% E_3 = load('data/eigen_mats/pbc_L_variable/values_81_3.mat').E(1:Neig);
% E_4 = load('data/eigen_mats/pbc_L_variable/values_81_4.mat').E(1:Neig);

x = [1:Neig];

% plot(x, E_0, 'x', x, 2*x*pi, x, E_1, 'x', x, 2*x*pi/(8/9), x, E_2, 'x', x, 2*x*pi/(64/81));


% figure;
%plot(x, E_1, x, E_2, x, E_3, x, E_4);
% plot(x, E_1, 'x', x, E_2, 'x', x, E_3, 'x', x, E_4, 'x');
% plot(x, E_1 - E_1(1), 'x', x/(3^2), E_2 - E_2(1), 'x', x/(9^2), E_3 - E_3(1), 'x', x/(27^2), E_4 - E_4(1), 'x');

% plot([1:100], E_2(1:100), 'x', x / 3^2 , E_3, 'x');
% plot(x, E_0, 'o', x/(3^2), E_1, 'x', x/(9^2), E_2, 'x', x/(27^2), E_3, 'x', x/(81^2), E_4, 'x');
% xlim([0, 40]);
% lgnd = legend('no potential', 'iteration 1', 'iteration 2', 'iteration 3', 'iteration 4'); 
% lgnd.Location = 'southeast'; 
% ylabel('Energy $\left ( \frac{\hbar ^2}{2m} \right )$','Interpreter','latex'); 
% xlabel('Eigenstate','Interpreter','latex'); 
% title('Energy spectrum for different iterations of Sierpinski carpet','Interpreter','latex'); 

% plot((1:20)*9, E_0, 'o', 1:20, E_1(1:20), 'x', (1:180)/(3^2), E_2(1:180), 'x', (1:1620)/(9^2), E_3(1:1620), 'x', x/(27^2), E_4, 'x');

% G_0 = gap_values(E_0);
% G_1 = gap_values(E_1);
% G_2 = gap_values(E_2);
% G_3 = gap_values(E_3);
% G_4 = gap_values(E_4);

% figure;
% plot([1:100], G_2(1:100) * 3^2, 'x', x / 3^2 , G_3 * 3^3, 'x');
% plot(x, G_1, 'x', x/3, G_2, 'x', x/(3^2), G_3, 'x', x/(3^3), G_4, 'x');


% -------------------------------------------------------------------------------------------
% FIGURE TFG REPORT 3281 eigenenergies 
% plot(x, E, 'x'); 
% plot(x, E_0, 'x', x, E_1, 'x', x, E_2,'x', x, E_3 ,'x', x, E_4,'x', x, (pi^2/2)*(4/pi)*x, x, ((pi^2/2)*(4/pi)*x)/(8/9)); 
% plot(x, E_0-E_0(1),'x', x/(3^2), E_1-E_1(1),'x', x/(9^2), E_2-E_2(1),'x', x/(27^2), E_3 -E_3(1),'x', x/(81^2), E_4-E_4(1),'x'); 

plot(x, E_0, 'x', x/(2^2), E_1, 'x', x/(4^2), E_2, 'x');
lgnd = legend('L = 1, N = 20', 'L = 2, N = 40', 'L = 4, N = 80'); 
lgnd.Location = 'southeast'; 
ylabel('Energy $\left ( \frac{\hbar ^2}{2m} \right )$','Interpreter','latex'); 
xlabel('Eigenstate','Interpreter','latex'); 
title('Eigenenergies without external potential, with x coordinate escalated','Interpreter','latex'); 
xlim([0, 25]);
% ylim([0, 500]);
% -------------------------------------------------------------------------------------------


function r = gap_values(E)
    V = zeros(length(E) - 1, 1);
    for j=(2:length(E))
        V(j) = E(j) - E(j - 1);
    end
    r = V;
end



% plot(x, E_0, 'o', x, E_1, 'x', x, E_2, '^', x, E_3, '+', x, E_4, '*');

% lgnd = legend('$\alpha = 0$', '$\alpha = 0.01$', '$\alpha = 0.05$', '$\alpha = 0.1$', '$\alpha = 0.2$', 'Interpreter','latex'); 
% lgnd.Location = 'southeast'; 
% ylabel('Energy $\left ( \frac{\hbar ^2}{2m} \right )$','Interpreter','latex'); 
% xlabel('Eigenstate','Interpreter','latex'); 
% title('Eigenenergies for fractal potential with noise','Interpreter','latex'); 
% xlim([0 40]);