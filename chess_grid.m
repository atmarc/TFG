function r = chess_grid(n, squares)
    % n -> size of the returning matrix
    % squares -> number of squares per row

    pot = 10000;
    M = zeros(n);

    for i=(1:n)
        for j=(1:n)
            x = mod(floor((i - 1)/(n/squares)), 2) == 0;
            y = mod(floor((j - 1)/(n/squares)), 2) == 0;
            if x && y || ~x && ~y
                M(i, j) = pot;
            end
        end
    end
    r = M;
end