function [S_valid] = testbimodality(S,vel,hartigan_check,hartigan_pthresh,Nboot)
%test_bimodality validate subset of vectors based on removing non-dominant modes in a multimodal distribution
%
%	Author: Abel-John Buchner
%
%   Input
%       S: logical array (1xN) flagging the subset of vectors to consider
%       vel: full array (3xN, only implemented 3D) of velocity (or position, or anything else really) vectors
%       hartigan_check: logical, 1= check unimodality with Hartigan's dip test, 0= assume bimodality without checking 
%       hartigan_pthresh: p-value threshold to apply to hartigan's dip statistic
%
%   Output
%       S_valid: logical array (1xN) flagging the vectors which passed the test
%
%   Note: Can only remove second mode in a bimodal distribution, cannot handle more modes

warning('testbimodal.m: TO BE CHECKED')

% correct input
if nargin<5
    Nboot = 250;
end
if nargin<4
    hartigan_pthresh = 0.05;
end
if nargin<3
    hartigan_check   = 1;
end

% find
S_indx = find(S);

% test Hartigan's dip statistic (Hartigan and Hartigan, 1985, Annals of Statistics)
vel_buf = vel(:,S);

% get dip statistic and p-value (p-value<0.05 implies 95% confidence ) that it is bimodal
if hartigan_check
    for ii = 1:3
        %[dipstat(ii), p_sigval(ii)] = hartigan_with_significance(vel_buf(ii,:));
        dipstat(ii) = HartigansDipTest(vel_buf(ii,:));
    end
    % calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
    % ie. we are testing for a variable that is random on some interval, how often it will produce a dip statistic higher than that yielded from the measured sample
    dipstat_boot = zeros(Nboot,1);
    vecsize = [1 size(vel_buf,2)];
    for ii = 1:Nboot % can I rewrite HartigansDipTest to take multiple vector inputs and so remove this loop?...
       unifpdf_boot     = sort(rand(vecsize));
       dipstat_boot(ii) = HartigansDipTest(unifpdf_boot);
    end
    p_sigval=sum(min(dipstat)<dipstat_boot)/Nboot;    
else
    p_sigval = 0; % set p-value to zero to trigger segmentation without hartigan's test
end

% if bimodal, only keep peak with larger population
if p_sigval <= hartigan_pthresh
    % separate data into 2 peaks (assuming only 2 for now). If doesn't converge, then keep all vectors
    try
        gmmod = fitgmdist(vel_buf',2,'MaxIter',75,'TolFun',1e-6); % works for any dimensionality data
        indxgrouping = cluster(gmmod,vel_buf')';
    catch
        indxgrouping = ones(1,size(vel_buf,2));
    end
    S_valid_indx = S_indx(indxgrouping==mode(indxgrouping)); % apply condition to all three components
else
    S_valid_indx = S_indx;
end

% write
S_valid = zeros(size(S));
S_valid(S_valid_indx) = 1;

end

