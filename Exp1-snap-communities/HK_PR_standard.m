% Need to add path to hkgrow, if not stored in 'algorithms' folder, then
% update the path here
addpath('../algorithms/hkgrow/')
addpath('../algorithms/pprpush')
addpath('datasets')

%%
datasets = {'DBLP','Amazon','LiveJournal','Orkut','Youtube'};
alpha = .99;
percentofset = 5;

% Tolerance parameters for PageRank Push
tols = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7];

for graph = datasets
    
    load(strcat(char(graph),'-top10.mat'))
    outputfile = strcat('Output',num2str(percentofset),'/',char(graph),'_hkpr_standard.txt');
    outputmatrix = strcat('Output',num2str(percentofset),'/',char(graph),'_hkpr_standard.mat');
    load(strcat(char(graph),'-seed-starter.mat'))

    if percentofset == 2
      S = S2;
    elseif percentofset == 3
      S = S3;
    else
      S = S5;
    end

    comm = C;
    d = sum(A,1)';
    volA = sum(nonzeros(A));
    n = size(A,1);
    numcom = size(comm,2);

    pr_sets = sparse(n,numcom);
    hk_sets = sparse(n,numcom);
    pr_stats = zeros(6,numcom);
    hk_stats = zeros(6,numcom);

    fid = fopen(outputfile,'w');
    for commID = 1:numcom

        % Target Community Stats
        Target = find(comm(:,commID));
        TarSize = numel(Target);
        [cutT,volT,edgesT,Tcond] = set_stats(A,Target,volA);
        Tsize = numel(Target);
        fprintf('\nCommunity %d has %d nodes and a conductance of %f \n',commID, Tsize,Tcond);
        fprintf(fid,'\nCommunity %d has %d nodes and a conductance of %f \n',commID, Tsize,Tcond);

        % Information about known target nodes
        R = find(S(:,commID));
        [cutR,volR,edgesR,condR] = set_stats(A,R,volA);
        [pR,rR,fR] = AdjustedPRF(Target,R,[]);
        fprintf('Seed set: conductance = %f, size = %d, pr = %f, re = %f, f1 = %f \n',condR,numel(R),pR,rR,fR);
        fprintf(fid,'Seed set: conductance = %f, size = %d, pr = %f, re = %f, f1 = %f \n',condR,numel(R),pR,rR,fR);

        fprintf('%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');
        fprintf(fid,'%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');

        % Use the full seed set, check Push output for many tolerance
        % parameters
        PR = [];
        condPR = 1;
        tic
        for tol = tols
            [aPR,acondPR,acutPR,avolPR,y] = pprpush_weighted_mex(A, d, R, tol, alpha,-1-1.e-4);
            if acondPR < condPR
                condPR = acondPR;
                PR = aPR;
            end
        end
        tPR = toc;

        [pPR,rPR,fPR] = AdjustedPRF(Target,PR,[]);
        fprintf('PR: \t %d \t %f \t %f \t %f \t %f \t %f \n',numel(PR),tPR,condPR,pPR,rPR,fPR)
        fprintf(fid,'PR: \t %d \t %f \t %f \t %f \t %f \t %f \n',numel(PR),tPR,condPR,pPR,rPR,fPR);

        pr_sets(PR,commID) = 1;
        if isnan(fPR)
         fPR = 0;
        end
        pr_stats(:,commID) = [tPR;numel(PR);condPR;pPR;rPR;fPR];

        % Do the same for HK-grow
        t_vals = [5 10 20 40 80];
        eps_vals = [1e-4 1e-4 1e-3 5*1e-3 1e-2];
        debugflag = 0;
        condHK = Inf;
        HK = [];

        tic
        for ei=1:numel(t_vals)
            [curset, condt,~,~,~,~] = hkgrow_mex(A, R, t_vals(ei), eps_vals(ei), debugflag);

            if condtemp < condHK
                condHK = condtemp;
                HK = curset;
                fHK = fHKt;
                pHK = pHKt;
                rHK = rHKt;
            end
        end
        tHK = toc;
        fprintf('HK: \t %d \t %f \t %f \t %f \t %f \t %f \n',numel(HK),tHK,condHK,pHK,rHK,fHK)
        fprintf(fid,'HK: \t %d \t %f \t %f \t %f \t %f \t %f \n',numel(HK),tHK,condHK,pHK,rHK,fHK);

        hk_sets(HK,commID) = 1;
        if isnan(fHK)
            fHK = 0;
        end
        hk_stats(:,commID) = [tHK;numel(HK);condHK;pHK;rHK;fHK];

    end
    
    save(outputmatrix,'pr_sets','pr_stats','hk_sets','hk_stats','tols','alpha')
end
