N = 3^9;

distance = zeros(1,100);
parfor i=1:9
    distance(i) = distance(i) + i;
end
