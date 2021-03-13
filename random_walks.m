global L;
global dt;
global Rmax;
global N;
global grid;

iterations = 5000;
n_walkers = 1000;
L = 1;
N = 3^6;
dt = L/(3^6); % L/3^i
Rmax = L/2;
 
grid = sierpinski(N, 6, true);
data = zeros(1,iterations);
walkers = init_walkers(n_walkers);

for i=1:n_walkers
    if mod(i, 50) == 0
        display(['Walker ' num2str(i)]);
    end
    
    for j=1:iterations
        walkers(i) = update_walker(walkers(i));
        dx = walkers(i).init_x - walkers(i).x;
        dy = walkers(i).init_y - walkers(i).y;

        data(j) = data(j) + dx^2 + dy^2;
        %display([num2str(g_posx) ' ' num2str(g_posy) ' ' num2str(walkers(i).x_pbc) ' ' num2str(walkers(i).y_pbc)]);
        
    end
end

%imagesc(grid);

% mean value
data = data / n_walkers;
x = 1:iterations;
x = x * dt;

p = polyfit(x, data, 1);
y_pred = polyval(p, x);

disp('Ajuste:');
disp(['y = x*' num2str(p(1)) ' + ' num2str(p(2))]);
plot(x, data, x, 2*x, x, y_pred);

function res = init_walkers(n_walkers)
    global L; global Rmax;
    walkers(n_walkers) = struct();
    for i=1:n_walkers
        x = rand() * L - Rmax;
        y = rand() * L - Rmax;
        
        while ~valid_pos(x, y)
            x = rand() * L - Rmax;
            y = rand() * L - Rmax;
        end    

        walkers(i).x = x;
        walkers(i).x_pbc = x;
        walkers(i).y = y;
        walkers(i).y_pbc = y;
        walkers(i).init_x = x;
        walkers(i).init_y = y;
    end
    res = walkers;
end

function res = update_walker(walker)
    tmp = get_random_pos(walker);
    
    while ~valid_pos(tmp.x_pbc, tmp.y_pbc)
        tmp = get_random_pos(walker);
    end

    res = tmp;
end

function r = valid_pos(x, y)
    global grid; global N; global L; global Rmax;
    
    g_posx = floor((x + Rmax) * N/L) + 1;
    g_posy = floor((y + Rmax) * N/L) + 1;
    r = grid(g_posx, g_posy) == 0;
end

function r = get_random_pos(walker)
    global dt; global Rmax; global L;
    
    % Update movement
    [a, b] = random_values();
    walker.x = walker.x + a * sqrt(dt);
    walker.x_pbc = walker.x_pbc + a * sqrt(dt);

    walker.y = walker.y + b * sqrt(dt);
    walker.y_pbc = walker.y_pbc + b * sqrt(dt);

    %disp([num2str(walker.x_pbc) ' ' num2str(walker.y_pbc)])
    % check if it is in a valid position
    if walker.x_pbc <= -Rmax
        walker.x_pbc = mod(walker.x_pbc, Rmax); 
    end
    if walker.x_pbc >= Rmax
        walker.x_pbc = mod(walker.x_pbc, -Rmax); 
    end
    if walker.y_pbc <= -Rmax 
        walker.y_pbc = mod(walker.x_pbc, Rmax); 
    end
    if walker.y_pbc >= Rmax
        walker.y_pbc = mod(walker.x_pbc, -Rmax);
    end
    %disp([num2str(walker.x_pbc) ' ' num2str(walker.y_pbc)])
    r = walker;
end

function [z1, z2] = random_values()
    % Box-Muller method
    phi = 2 * pi * rand();
    r = sqrt(-2 * log(rand()));
    
    z1 = r * cos(phi);
    z2 = r * sin(phi);
end