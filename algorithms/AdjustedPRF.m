function [pr,re,F1] = AdjustedPRF(True, Returned, Seed)
%
% Find the F1 detection score for a certain method, but only measured on
% the data that wasn't a part of the seed set.
%

Returned = setdiff(Returned,Seed);
Target = setdiff(True,Seed);

TruePos = intersect(Returned,Target);
pr = numel(TruePos)/numel(Returned);
re = numel(TruePos)/numel(Target);
F1 = 2*(pr*re)/(pr+re);
