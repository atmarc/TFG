A = csvread('track_mat.txt');
imagesc(A);
pause;

rec_lvl = 5;
iterations = 10^8;
data_size = 10^4;
L = 1;
distances = csvread('distances.txt');
dt = ((L/3^rec_lvl)/10)^2;

x = 1:(data_size + 1);
x = x * sqrt(dt) * iterations/data_size;

p = polyfit(x, distances, 1);
y_pred = polyval(p, x);

disp(['y = x*' num2str(p(1)) ' + ' num2str(p(2))]);

plot(x, distances, x, y_pred);