addpath('Output5')
% clear
% clc
datasets = {'DBLP','Amazon','LiveJournal','Orkut'};

for graph = datasets
   
  load(strcat(char(graph),'_crd.mat'))
  load(strcat(char(graph),'_hkpr_standard.mat'))
  load(strcat(char(graph),'slfs_0.1epsnosotf.mat'))
  
  % Stats format:
  % [time size conductance precision recall F1]';
  
  % Conductance averages
  crd_cond = mean(crd_stats(3,:));
  hk_cond = mean(hk_stats(3,:));
  pr_cond = mean(pr_stats(3,:));
  sl_cond = mean(sl_stats(3,:));
  fs_cond = mean(fs_stats(3,:));
  
  % F1 averages
  crd_f1 = mean(crd_stats(6,:));
  hk_f1 = mean(hk_stats(6,:));
  pr_f1 = mean(pr_stats(6,:));
  sl_f1 = mean(sl_stats(6,:));
  fs_f1 = mean(fs_stats(6,:));
  
  % precision averages
  crd_pr = mean(crd_stats(4,:));
  hk_pr = mean(hk_stats(4,:));
  pr_pr = mean(pr_stats(4,:));
  sl_pr = mean(sl_stats(4,:));
  fs_pr = mean(fs_stats(4,:));
  
  % recall averages
  crd_re = mean(crd_stats(5,:));
  hk_re = mean(hk_stats(5,:));
  pr_re = mean(pr_stats(5,:));
  sl_re = mean(sl_stats(5,:));
  fs_re = mean(fs_stats(5,:));
  
  % Time averages
  crd_time = mean(crd_stats(1,:));
  hk_time = mean(hk_stats(1,:));
  pr_time = mean(pr_stats(1,:));
  sl_time = mean(sl_stats(1,:));
  fs_time = mean(fs_stats(1,:));
  
  % size averages
  crd_size = round(mean(crd_stats(2,:)));
  hk_size = round(mean(hk_stats(2,:)));
  pr_size = round(mean(pr_stats(2,:)));
  sl_size = round(mean(sl_stats(2,:)));
  fs_size = round(mean(fs_stats(2,:)));
  
  fprintf(' \\midrule \n')
  fprintf('%s & \\alg{HK-relax} & %d & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n',char(graph),hk_size,hk_cond,hk_time,hk_pr,hk_re,hk_f1)
  fprintf(' & \\alg{Push} & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',pr_size,pr_cond,pr_time,pr_pr,pr_re,pr_f1)
  fprintf(' & \\alg{CRD}  & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',crd_size,crd_cond,crd_time,crd_pr,crd_re,crd_f1)
  fprintf(' & \\alg{SimpleLocal} & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',sl_size,sl_cond,sl_time,sl_pr,sl_re,sl_f1)
  fprintf(' & \\alg{FlowSeed} & %d & %.3f & %.3f & %.3f& %.3f & %.3f \\\\\n',fs_size,fs_cond,fs_time,fs_pr,fs_re,fs_f1)

  % Uncomment below for shorter version of the table
%   fprintf('%s & $\\phi$ & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',char(graph), fs_cond, sl_cond, hk_cond, pr_cond, crd_cond)
%   fprintf('& size & %d & %d & %d & %d & %d\\\\\n', fs_size, sl_size, hk_size, pr_size, crd_size)
%   fprintf('& F1 & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', fs_f1, sl_f1, hk_f1, pr_f1, crd_f1)
%   fprintf('& run. & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', fs_time, sl_time, hk_time, pr_time, crd_time)

  
end
  