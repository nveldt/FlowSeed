function [pr,re,F1] = PRF(Target, Returned)


TruePos = intersect(Returned,Target);
pr = numel(TruePos)/numel(Returned);
re = numel(TruePos)/numel(Target);
F1 = 2*(pr*re)/(pr+re);
