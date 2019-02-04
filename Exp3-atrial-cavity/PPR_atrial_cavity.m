% data/MRI_Graph contains: 'mri_A','mri_vol', 'cavity'

load('data/MRI_Graph.mat')

G = mri_A;
volA = mri_vol;
Target = cavity;
d = full(sum(G,2));
n = size(G,1);

load('data/cavity_seeds') % contains 'seed_start', 100 nodes from cavity

%% Run PPR push algorithm

addpath('../algorithms/pprpush')
alpha = .6;
tols = [1e-8, 1e-9, 1e-10, 1e-11,1e-12, 1e-13, 1e-14];

outputfile = strcat('PPR_Output/alpha_',num2str(alpha*100),'.txt');
fid = fopen(outputfile,'w');
        
sets = sparse(n,7);
stats= zeros(7,7);

% The full seed set consists of the 100 start nodes, plus neighbors
R = neighborhood(G,seed_start,2);
[cut,volR,edges,condR] = set_stats(G,R,volA);

% Report statistics for the seed set
fprintf('\tSeed set:conductance = %f, size = %d \n',condR,numel(R))
fprintf('%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');

fprintf(fid,'\tSeed set:conductance = %f, size = %d \n',condR,numel(R));
fprintf(fid,'%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');

% Run PPR push for a range of different tolerance values
for t = 1:numel(tols)
    tol = tols(t);
    tic
    [PR,cond,cut,vol,y] = pprpush_weighted_mex(G, d, R, tol, alpha,-1-1.e-4);
    tPR = toc;
    [pPR,rPR,fPR] = PRF(Target,PR);
    fprintf('%d \t %d \t %f \t %f \t %f \t %f \t %f \n',round(log10(tol)),numel(PR),tPR,cond,pPR,rPR,fPR)
    fprintf(fid,'%d \t %d \t %f \t %f \t %f \t %f \t %f \n',round(log10(tol)),numel(PR),tPR,cond,pPR,rPR,fPR);

    % save the output
    sets(PR,t) = 1;
    stats(:,t) = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];
end
fclose(fid);
save(strcat('PPR_Output/alpha_',num2str(alpha*100),'sets_stats'),'sets','stats')
