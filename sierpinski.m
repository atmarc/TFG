
function r = sierpinski(n, d)
  % Function to generate a sierpinski carpet of size nxn with d recursion levels
  M = zeros(n,n);
  display(['Generating sierpinski carpet of size: ' num2str(n) 'x' num2str(n)]);
  for i = 1:n
  % ------ To see execution -----
    if mod(i, 50) == 0
      display(['Rows computed: ' num2str(i) '/' num2str(n)]);
    end
  % -----------------------------
 
    for j = 1:n
      M(i,j) = sierpinski_rec(i - 1,j - 1, n, n, 0, d);
    end
  end
  r = M;
end


function v = sierpinski_rec(x, y, width, height, d, max_depth)
  % If is in the center square
  x2 = floor(x * 3 / width);
  y2 = floor(y * 3 / height);
  
  if (x2 == 1 && y2 == 1)
    v = 10000; % full
    
  else
    % Recursive case
    if d < max_depth
      x -= floor((x2 * width + 2) / 3);
      y -= floor((y2 * height + 2) / 3);
      width = floor((width + 2 - x2) / 3);
      height = floor((height + 2 - y2) / 3);
      v = sierpinski_rec(x, y, width, height, d + 1, max_depth);
    
    else
      v = 0; % empty
      
    end
  end  
end