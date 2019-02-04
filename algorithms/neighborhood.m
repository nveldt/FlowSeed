function N = neighborhood(A,R,k)
% N = neighborhood(A,R,k)
% Return the k neighborhood of seed set R in graph A
% This code is inefficient, but it isn't used too often.

n = size(A,1);
eR = zeros(n,1);
eR(R) = 1;
eN = eR;

for i = 1:k
    eN = A*eN > 0;
end

N = find(eN);

end