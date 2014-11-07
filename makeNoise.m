function NOISE = makeNoise(SIZE, LIM, CONF)

% Uncertainty bounds
Bounds = norminv( [(1 - CONF)/2, CONF + ((1 - CONF)/2)], 0, 1);

% Std dev
stdev = LIM / max(Bounds);

% Noise matrix
NOISE = stdev * randn(SIZE);

end

