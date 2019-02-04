datasets = {'DBLP','Amazon','LiveJournal','Orkut'};

clc
for graph = datasets
   
  load(strcat(char(graph),'_crd.mat'))
  
  load(strcat(char(graph),'_hkpr.mat'))
  hk_sets_1 = hk_sets;
  hk_stats_1 = hk_stats;
  
  load(strcat(char(graph),'_hkpr_range.mat'))
  pr_sets_1 = pr_sets;
  pr_stats_1 = pr_stats;

  load(strcat(char(graph),'slfs_0.1epsnosotf.mat'))
  fs_sets_10 = fs_sets;
  fs_stats_10 = fs_stats;
  sl_sets_10 = sl_sets;
  sl_stats_10 = sl_stats;
  
  % Stats format:
  % [time size conductance precision recall F1]';
  
  % Conductance averages
  crd_cond = mean(crd_stats(3,:));
  hk_cond = mean(hk_stats_1(3,:));
  pr_cond = mean(pr_stats_1(3,:));
  sl_cond = mean(sl_stats_10(3,:));
  fs_cond = mean(fs_stats_10(3,:));
  
  % F1 averages
  crd_f1 = mean(crd_stats(6,:));
  hk_f1 = mean(hk_stats_1(6,:));
  pr_f1 = mean(pr_stats_1(6,:));
  sl_f1 = mean(sl_stats_10(6,:));
  fs_f1 = mean(fs_stats_10(6,:));
  
  % precision averages
  crd_pr = mean(crd_stats(4,:));
  hk_pr = mean(hk_stats_1(4,:));
  pr_pr = mean(pr_stats_1(4,:));
  sl_pr = mean(sl_stats_10(4,:));
  fs_pr = mean(fs_stats_10(4,:));
  
  % recall averages
  crd_re = mean(crd_stats(5,:));
  hk_re = mean(hk_stats_1(5,:));
  pr_re = mean(pr_stats_1(5,:));
  sl_re = mean(sl_stats_10(5,:));
  fs_re = mean(fs_stats_10(5,:));
  
  % Time averages
  crd_time = mean(crd_stats(1,:));
  hk_time = mean(hk_stats_1(1,:));
  pr_time = mean(pr_stats_1(1,:));
  sl_time = mean(sl_stats_10(1,:));
  fs_time = mean(fs_stats_10(1,:));
  
  % size averages
  crd_size = round(mean(crd_stats(2,:)));
  hk_size = round(mean(hk_stats_1(2,:)));
  pr_size = round(mean(pr_stats_1(2,:)));
  sl_size = round(mean(sl_stats_10(2,:)));
  fs_size = round(mean(fs_stats_10(2,:)));
  
  fprintf(' \\midrule \n')
  fprintf('%s & \\alg{HK-relax} & %d & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n',char(graph),hk_size,hk_cond,hk_time,hk_pr,hk_re,hk_f1)
  fprintf(' & \\alg{Push} & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',pr_size,pr_cond,pr_time,pr_pr,pr_re,pr_f1)
  fprintf(' & \\alg{CRD}  & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',crd_size,crd_cond,crd_time,crd_pr,crd_re,crd_f1)
  fprintf(' & \\alg{SimpleLocal} & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',sl_size,sl_cond,sl_time,sl_pr,sl_re,sl_f1)
  fprintf(' & \\alg{FlowSeed} & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',fs_size,fs_cond,fs_time,fs_pr,fs_re,fs_f1)

  
end
  