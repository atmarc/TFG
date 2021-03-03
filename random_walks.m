L = 1;
Rmax = L/2;
N = 300;
dt = L/N;
iterations = 5000;
n_walkers = 100;

props.N = N; props.L = L; props.dt = dt; props.Rmax = Rmax;

grid = zeros(N,N);
data = zeros(n_walkers, iterations);

walkers = init_walkers(n_walkers);

for i=1:n_walkers
    for j=1:iterations
        walkers(i) = update_walker(walkers(i), props);
        w = walkers(i);
        
        %display([num2str(w.g_posx) ' ' num2str(w.g_posy)]);
        grid(w.g_posx, w.g_posy) = grid(w.g_posx, w.g_posy) + 1;
        data(i, j) = w.x_pbc^2 + w.y_pbc^2;
    end
end

%imagesc(grid);
x = 1:iterations;
y = zeros(1,iterations);
for i=1:iterations
    y(i) = mean(data(:, i));
end

plot(y);
%plot(x,data(1,:),x,data(2,:),x,data(3,:),x,data(4,:),x,data(5,:));




function res = init_walkers(n_walkers)
    walkers(n_walkers) = struct();
    for i=1:n_walkers
        walkers(i).x = 0;
        walkers(i).x_pbc = 0;
        walkers(i).y = 0;
        walkers(i).y_pbc = 0;
        walkers(i).g_posx = -1;
        walkers(i).g_posy = -1;
    end
    res = walkers;
end

function res = update_walker(walker, props)
    dt = props.dt; N = props.N; L = props.L; Rmax = props.Rmax;
    
    [a, b] = random_values();
    walker.x = walker.x + a * dt;
    walker.x_pbc = walker.x_pbc + a * dt;
    walker.y = walker.y + b * dt;
    walker.y_pbc = walker.y_pbc + a * dt;
    
    walker.g_posx = round((walker.x + Rmax) * N/L) + 1;
    walker.g_posy = round((walker.y + Rmax) * N/L) + 1;
    
    % not valid position
    if walker.g_posx <= 0 
        walker.x = walker.x + L;
        walker.g_posx = walker.g_posx + N; 
    end
    if walker.g_posx > N 
        walker.x = walker.x - L;
        walker.g_posx = walker.g_posx - N; 
    end
    if walker.g_posy <= 0 
        walker.y = walker.y + L;
        walker.g_posy = walker.g_posy + N; 
    end
    if walker.g_posy > N 
        walker.y = walker.y - L;
        walker.g_posy = walker.g_posy - N; 
    end
    res = walker;
end

function [z1, z2] = random_values()
    % Box-Muller method
    phi = 2 * pi * rand();
    r = sqrt(-2 * log(rand()));
    
    z1 = r * cos(phi);
    z2 = r * sin(phi);
end