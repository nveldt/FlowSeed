%%
%A = load_graph('four-clusters');

addpath('~/data/MatlabGraphs/')
load lesmis
A = Problem.A;
n = size(A,1);
d = full(sum(A,2));
alpha = 0.85;
tol = 1e-4;
v = zeros(n,1);
v(1) = 1;

P = diag(sparse(1./d))*A;

%%
[xr,rr] = acl_method(P, (1-alpha)*v, d, alpha, tol, 1);

%% Fix bad termination criteria
[xr,rr] = acl_method_1(P, (1-alpha)*v, d, alpha, tol);

%%
% [bestset,cond,cut,vol,y] = pprpush_weighted_mex(A,degs,set,eps,alpha)
[bestset,cond,cut,vol,y] = pprpush_weighted_mex(A, d, 1, tol, alpha);

%%
fprintf('Simple test\n');
norm(y-xr)
assert(norm(y-xr) < 1e-15, 'simple test failed');

%% More exhaustive test
for i=1:n
    vv = zeros(n,1);
    vv(i) = 1;
    [xr,rr] = acl_method_1(P, (1-alpha)*vv, d, alpha, tol, 1); 
    [bestset,cond,cut,vol,y] = pprpush_weighted_mex(A, d, i, tol, alpha);
    if norm(y-xr) > n*eps(1)
        fprintf('Failed on four-clusters node %i\n', i);
        assert(false);
    end
end

%% Non-uniform weight test
%%
A = load_graph('four-clusters');
A = sprand(A);
A = triu(A,1);
A = A + A';
n = size(A,1);
d = full(sum(A,2));
alpha = 0.85;
tol = 1e-4;
v = zeros(n,1);
v(1) = 1;

P = diag(sparse(1./d))*A;

%%
[xr,rr] = acl_method(P, (1-alpha)*v, d, alpha, tol, 1);

%%
[xr,rr] = acl_method_1(P, (1-alpha)*v, d, alpha, tol);

%%
% [bestset,cond,cut,vol,y] = pprpush_weighted_mex(A,degs,set,eps,alpha)
[bestset,cond,cut,vol,y] = pprpush_weighted_mex(A, d, 1, tol, alpha);

%%
fprintf('Non-uniform weight test\n');
norm(y-xr)
assert(norm(y-xr) < 1e-15, 'non-uniform weighted test failed');

%% Larger test
n = 1000;
A = sprand(n,n,10/1000);
A = triu(A,1);
A = A + A';
d = full(sum(A,2));
P = diag(sparse(1./d))*A;
alpha = 0.85;
tol = 1e-3;
v = zeros(n,1);
v(1) = 1;
%%
%[xr,rr] = acl_method_1(P, (1-alpha)*v, d, alpha, tol, 1); 
[bestset,cond,cut,vol,y] = pprpush_weighted_mex(A, d, 1, tol, alpha);
fprintf('Simple large test\n');
norm(y-xr)
assert(norm(y-xr) < 1e-15, 'simple large test failed');

%%
for i=1:n
    vv = zeros(n,1);
    vv(i) = 1;
    [xr,rr] = acl_method_1(P, (1-alpha)*vv, d, alpha, tol, 1); 
    [bestset,cond,cut,vol,y] = pprpush_weighted_mex(A, d, i, tol, alpha);
    if norm(y-xr) > n*eps(1)
        fprintf('Failed on random graph test with node %i\n', i);
        assert(false);
    end
end

%% Now compare with multiple seeds

n = 1000;
A = sprand(n,n,10/1000);
A = triu(A,1);
A = A + A';
d = full(sum(A,2));
P = diag(sparse(1./d))*A;
alpha = 0.85;
tol = 1e-3;
v = zeros(n,1);
v(1) = 1/2;
v(32) = 1/2;

[xr,rr] = acl_method_1(P, (1-alpha)*v, d, alpha, tol, 1); 
[bestset,cond,cut,vol,y] = pprpush_weighted_mex(A, d, sparse(v), tol, alpha);
fprintf('Multiple seed test\n');
norm(y-xr)

%%
for i=1:n
    nz = randi(100,1);
    vv = sprand(n,1,nz/n);
    vv = vv/sum(vv);
    [xr,rr] = acl_method_1(P, (1-alpha)*vv, d, alpha, tol, 1); 
    fprintf('Test %i\n', i);
    [bestset,cond,cut,vol,y] = pprpush_weighted_mex(A, d, vv, tol, alpha);
    if norm(y-xr) > n*eps(1)
        fprintf('Failed on node %i\n', i);
        assert(false);
    end
end
fprintf('Passed multiple seed test\n');

%% Test prl1_push vs. prl1_cvx
A = load_graph('four-clusters');
A = sprand(A);
A = triu(A,1);
A = A + A';
n = size(A,1);
d = full(sum(A,2));
alpha = 0.85;
tol = 1e-4;
v = zeros(n,1);
v(1) = 1;

for i=1:n
    vv = zeros(n,1);
    vv(i) = 1;
    cutalpha = (1-alpha)/alpha;
    cutkappa = d(i)*tol/alpha;
    fprintf('prl1-cvx test %i\n', i);
    zpush = prl1_push(A, cutalpha, vv, d, cutkappa);
    zcvx = prl1_cvx(A, cutalpha, vv, d, cutkappa);
    if norm(zpush-zcvx,inf) > 1e-4 % CVX solves to a tolerance of 1e-8
        fprintf('Failed on node %i\n', i);
        assert(false);
    end
end
fprintf('Passed prl1 vs. cvx test\n');


%% 
% Multiple seeds
cvx_precision high
for i=1:n
    nz = randi(15,1);
    vv = spones(sprand(n,1,nz/n));
    cutalpha = (1-alpha)/alpha;
    cutkappa = d(i)*tol/alpha;
    zpush = prl1_push(A, cutalpha, vv, d, cutkappa);
    zcvx = prl1_cvx(A, cutalpha, vv, d, cutkappa);
    if norm(zpush-zcvx,inf) > 2e-6 % CVX solves to a tolerance of 1e-6 or so
        fprintf('Failed on node %i\n', i);
        assert(false);
    end
end
cvx_precision([]);

%% 
% Multiple seeds with uneven weights
cvx_precision high
for i=1:n
    nz = randi(15,1);
    vv = sprand(n,1,nz/n);
    cutalpha = (1-alpha)/alpha;
    cutkappa = d(i)*tol/alpha;
    zpush = prl1_push(A, cutalpha, vv, d, cutkappa);    
    zcvx = prl1_cvx(A, cutalpha, vv, d, cutkappa);
    if norm(zpush-zcvx,inf) > 1e-4 % CVX solves to a tolerance of 1e-8
        fprintf('Failed on node %i\n', i);
        assert(false);
    end
end
cvx_precision([]);

%% Compare against gurobi
fprintf('If you don''t have gurobi, these tests will fail.\n');

%% Compare against prl1_gurobi
A = load_graph('four-clusters');
A = sprand(A);
A = triu(A,1);
A = A + A';
n = size(A,1);
d = full(sum(A,2));
alpha = 0.85;
tol = 1e-4;
v = zeros(n,1);
v(1) = 1;

for i=1:n
    vv = zeros(n,1);
    vv(i) = 1;
    cutalpha = (1-alpha)/alpha;
    cutkappa = d(i)*tol/alpha;
    zpush = prl1_push(A, cutalpha, vv, d, cutkappa);
    zgurobi = prl1_gurobi(A, cutalpha, vv, d, cutkappa);
    if norm(zpush-zgurobi,inf) > 1e-4 % Gurobi solves to a tolerance of 1e-8
        fprintf('Failed on node %i\n', i);
        assert(false);
    end
end

%%

for i=1:n
    vv = zeros(n,1);
    vv(i) = 1;
    cutalpha = (1-alpha)/alpha;
    cutkappa = d(i)*tol/alpha;
    zpush = prl1_push(A, cutalpha, vv, d, cutkappa);
    zgurobi = prl1_gurobi(A, cutalpha, vv, d, cutkappa);
    if norm(zpush-zgurobi,inf) > 1e-4 % Gurobi solves to a tolerance of 1e-8
        fprintf('Failed on node %i\n', i);
        assert(false);
    end
end

%%
% Larger test
n = 1000;
A = sprand(n,n,10/1000);
A = triu(A,1);
A = A + A';
d = full(sum(A,2));
P = diag(sparse(1./d))*A;
alpha = 0.85;
tol = 1e-3;
tpush = 0;
tgurobi = 0;
for i=1:n
    vv = zeros(n,1);
    vv(i) = 1;
    cutalpha = (1-alpha)/alpha;
    cutkappa = d(i)*tol/alpha;
    t0 = tic;
    zpush = prl1_push(A, cutalpha, vv, d, cutkappa);
    tpush = tpush + toc(t0);
    t0 = tic;
    zgurobi = prl1_gurobi(A, cutalpha, vv, d, cutkappa);
    tgurobi = tgurobi + toc(t0);
    if norm(zpush-zgurobi,inf) > 1e-4 % Gurobi solves to a tolerance of 1e-8
        fprintf('Failed on node %i\n', i);
        assert(false);
    end
    if mod(i,10) == 1, fprintf('node %4i tpush = %8.1f tgurobi = %8.1f\n', i, tpush, tgurobi); end
end
