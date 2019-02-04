rfunction x = pprpush(A,d,beta,v,tol)
% PPRPUSH Solve a personalized PageRank problem using the ACL Push method
% 

n = size(A,1);
if numel(v) == n
    v = sparse(v);
    if any ( v > 1 ) 
        assert( false, 'not yet implemented, please specify v as a vector');
    end
    sumv = sum(v);
    assert ( abs(sumv - 1) < 10*nnz(v)*eps(1), 'v is not stochastic');
else
    % v must be a set 
end

[~,~,~,~,x] = pprpush_weighted_mex(A, d, v, tol, beta);
