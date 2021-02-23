
function r = sierpinski(n, d, verbose, potencial)
  % Function to generate a sierpinski carpet of size nxn with d recursion levels
  if nargin < 3
    verbose = false;  
  end
  
  if nargin < 4
    potencial = 10000;  
  end
  
  M = zeros(n,n);
  
  if verbose
    display(['Generating sierpinski carpet of size: ' num2str(n) 'x' num2str(n)]);
  end
  
  for i = 1:n
  % ---------- Verbose ----------
    if (verbose && mod(i, 50) == 0)
      display(['Rows computed: ' num2str(i) '/' num2str(n)]);
    end
  % -----------------------------
 
    for j = 1:n
      M(i,j) = sierpinski_rec(i - 1,j - 1, n, n, 1, d, potencial);
    end
  end
  r = M;
end


function v = sierpinski_rec(x, y, width, height, d, max_depth, potencial)
  % If is in the center square
  x2 = floor(x * 3 / width);
  y2 = floor(y * 3 / height);
  
  if (x2 == 1 && y2 == 1)
    v = potencial; % full
    
  else
    % Recursive case
    if d < max_depth
      x = x - floor((x2 * width + 2) / 3);
      y = y - floor((y2 * height + 2) / 3);
      width = floor((width + 2 - x2) / 3);
      height = floor((height + 2 - y2) / 3);
      v = sierpinski_rec(x, y, width, height, d + 1, max_depth, potencial);
    
    else
      v = 0; % empty
      
    end
  end  
end