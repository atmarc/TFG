thresholds = 0.01:0.05:3;
[n_values, ~] = size(thresholds);
dimensions = zeros(n_values, 1);
index = 1;
for x=thresholds
    area = PSI_2 > x;
    Mat = area(1:364, 1:364);
    dim = b_fractal_dimension(Mat);
    %disp(dim);
    dimensions(index) = dim;
    index = 1 + index;
end

%figure;
plot(thresholds, dimensions, 'x'); xlabel("Wavefunction threshold"); ylabel("Box-counting dimension");
title("Box-counting dimension of the wavefunction depending on a threshold");
function dim = b_fractal_dimension(Mat)
    verbose = false;

    data = box_counting_method(Mat);
    Y = log(data(:, 1));
    X = log(data(:, 2));

    p = polyfit(X, Y, 1);
    x_pred = [min(X) max(X)];
    y_pred = polyval(p, x_pred);

    if verbose
        disp(['Dimension: ' num2str(p(2))]);
        plot(X, Y, 'x', x_pred, y_pred); xlabel("log(N)"); ylabel("log(s)"); title("Box-counting dimension of the fractal");
    end
    
    dim = p(2);
end


function r = box_counting_method(Mat)
    [rows, ~] = size(Mat);
    scale = 2:40;
    [~, n_data] = size(scale);
    data = zeros(n_data, 2);
    
    for j=1:n_data
        e = floor(rows/scale(j));
        [m, n] = size(Mat);
        N = 0;
        for x=1:e:m
            for y=1:e:n
                if is_border(Mat, x, y, e)
                    N = N + 1;
                end
            end        
        end
        data(j, 1) = N;
        data(j, 2) = j;
        %disp([num2str(N) ' ' num2str(scale(j))]);
    end
    r = data;
end


function r = is_border(Mat, fila, col, e)
    empty = false;
    fill = false;
    [m, n] = size(Mat);
%     global prova;
    
    for x=0:(e - 1)
        for y=0:(e - 1)
            if fila + x <= m && col + y <= n && Mat(fila + x, col + y) == 0
                empty = true;
            else
                fill = true;
            end
        end
    end
   
    % ----------------------------
%     if empty && fill
%         for x=0:(e - 1)
%             for y=0:(e - 1)
%                 if fila + x <= m && col + y <= n
%                     prova(fila + x, col + y) = 1;
%                 end
%             end
%         end
%     end
    % ----------------------------


    r = empty && fill;
end