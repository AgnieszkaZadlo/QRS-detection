function signal_out = dilatation(signal,element)

K = length(signal);
M = length(element);
P = M - 1;

signal_prim = zeros(1, K+P);
signal_prim(M:end) = signal;
signal_prim(1:(M-1)) = repmat(signal(1), 1, M-1);

signal_out = zeros(1, K+P);

for k=M:K+P
    vect = zeros(1,M);
    for m=1:M
        vect(m) = signal_prim(k-m+1) + element(m);
    end
    signal_out(k) = max(vect);
end

signal_out=signal_out(M:end);

end