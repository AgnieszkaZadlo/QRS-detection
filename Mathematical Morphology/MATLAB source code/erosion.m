function signal_out = erosion(signal,element)

K = length(signal);
M = length(element);
P = M - 1;

signal_prim = zeros(1, K+P);
signal_prim(1:K) = signal;
signal_prim(K+1:end) = repmat(signal(end), 1, M-1);

signal_out = zeros(1, K+P);

for k=1:K
    vect = zeros(1,M);
    for m=1:M
        vect(m) = signal_prim(k+m-1) - element(m);
    end
    signal_out(k) = min(vect);
end

signal_out=signal_out(1:K);

end