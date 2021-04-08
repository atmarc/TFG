global L; global dt; global Rmax; global N; global grid;
% So random numbers do not repeat
rng('shuffle');

iterations = 10000;
n_walkers = 2000;
rec_lvl = 1;
L = 1;
Rmax = L/2;
N = 3^9;
data_range = iterations;


for i_rec=1:9
    rec_lvl = i_rec;
    disp(['REC LVL: ' num2str(i_rec)]);
    
    % Load fractal
    grid = matfile(['fractals/sierpinski_' num2str(i_rec) '_' num2str(N) '.mat']).grid;
    disp('sierpinski loaded');
    
    dt = (L/(3^rec_lvl))^2/10; % L/3^i
    distances = zeros(1, data_range);
    walkers = init_walkers(n_walkers);

    
    for i=1:n_walkers
        if mod(i, 100) == 0
            display(['Walker ' num2str(i)]);
        end
        
        i_data = 1;
        for j=1:iterations
            if mod(j, 50000) == 0
                display(['Iteration ' num2str(j)]);
            end
            walkers(i) = update_walker(walkers(i));
            dx = walkers(i).init_x - walkers(i).x;
            dy = walkers(i).init_y - walkers(i).y;

            if mod(j, floor(iterations/data_range)) == 0 && i_data <= data_range
                distances(i_data) = distances(i_data) + dx^2 + dy^2;
                i_data = i_data + 1;
            end
        end
    end


    % ------------------ PLOTS ------------------
    % mean value
    distances = distances / n_walkers;

    x = 1:data_range;
    x = x * dt * iterations/data_range;

    p = polyfit(x, distances, 1);
    y_pred = polyval(p, x);

    disp('Ajuste:');
    disp(['y = x*' num2str(p(1)) ' + ' num2str(p(2))]);

    % plot(x, distances, x, 2*x, x, y_pred);
    % pause;
    save_to_file('random_walks_data.txt', [num2str(i_rec) ' ' num2str(p(1)) ' ' num2str(p(2))]);
end

% ----------------- METHODS -----------------

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
    
    if valid_pos(tmp.x_pbc, tmp.y_pbc)
        res = tmp;
    else
        res = walker;
    end
end

function r = valid_pos(x, y)
    global grid; global N; global L; global Rmax;
    
    g_posx = floor((x + Rmax) * N/L) + 1;
    g_posy = floor((y + Rmax) * N/L) + 1;
    r = grid(g_posx, g_posy) == 0;
end

function r = get_random_pos(walker)
    global dt; global Rmax;
    
    % Update movement
    [a, b] = random_values();
    walker.x = walker.x + a * sqrt(dt);
    walker.x_pbc = walker.x_pbc + a * sqrt(dt);

    walker.y = walker.y + b * sqrt(dt);
    walker.y_pbc = walker.y_pbc + b * sqrt(dt);
    
    % check if it is in a valid position to preserve pbc
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
    
    r = walker;
end

function [z1, z2] = random_values()
    % Box-Muller method
    phi = 2 * pi * rand();
    r = sqrt(-2 * log(rand()));
    
    z1 = r * cos(phi);
    z2 = r * sin(phi);
end
