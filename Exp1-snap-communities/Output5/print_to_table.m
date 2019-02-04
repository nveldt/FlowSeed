addpath('Output5')
datasets = {'DBLP','Amazon','LiveJournal','Orkut'};

clc
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
  
  fprintf('%s & $\\phi$ & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',char(graph), fs_cond, sl_cond, hk_cond, pr_cond, crd_cond)
  fprintf('& size & %d & %d & %d & %d & %d\\\\\n', fs_size, sl_size, hk_size, pr_size, crd_size)
  fprintf('& F1 & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', fs_f1, sl_f1, hk_f1, pr_f1, crd_f1)
  fprintf('& run. & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', fs_time, sl_time, hk_time, pr_time, crd_time)

  fprintf('\\midrule \n ')
  
end
  