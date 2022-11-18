function d = computeCohen_d(x1, x2, varargin)
% 
% call: d = computeCohen_d(x1, x2, varargin)
% 
% EFFECT SIZE of the difference between the two 
% means of two samples, x1 and x2 (that are vectors), 
% computed as "Cohen's d". 
% 
% If x1 and x2 can be either two independent or paired 
% samples, and should be treated accordingly:
%  
%   d = computeCohen_d(x1, x2, 'independent');  [default]
%   d = computeCohen_d(x1, x2, 'paired');
% 
% Note: according to Cohen and Sawilowsky:
%
%      d = 0.01  --> very small effect size
%      d = 0.20  --> small effect size
%      d = 0.50  --> medium effect size
%      d = 0.80  --> large effect size
%      d = 1.20  --> very large effect size
%      d = 2.00  --> huge effect size
%
%
% Ruggero G. Bettinardi (RGB)
% Cellular & System Neurobiology, CRG
% -------------------------------------------------------------------------------------------
%
% Copyright (c) 2017, Ruggero G. Bettinardi
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% Code History:
%
% 25 Jan 2017, RGB: Function is created
  
if nargin < 3, testType = 'independent'; 
else           testType = varargin{1}; 
end

% basic quantities:
n1       = numel(x1);
n2       = numel(x2);
mean_x1  = nanmean(x1);
mean_x2  = nanmean(x2);
var_x1   = nanvar(x1);
var_x2   = nanvar(x2);
meanDiff = (mean_x1 - mean_x2);

% select type of test:
isIndependent = strcmp(testType, 'independent');
isPaired      = strcmp(testType, 'paired');

% compute 'd' accordingly:
if isIndependent
    
    sv1      = ((n1-1)*var_x1);
    sv2      = ((n2-1)*var_x2);
    numer    =  sv1 + sv2;
    denom    = (n1 + n2 - 2);
    pooledSD =  sqrt(numer / denom); % pooled Standard Deviation
    s        = pooledSD;             % re-name
    d        =  meanDiff / s;        % Cohen's d (for independent samples)
    
elseif isPaired
    
    haveNotSameLength = ~isequal( numel(x1), numel(x2) );
    if haveNotSameLength, error('In a paired test, x1 and x2 have to be of same length!'), end
    
    deltas   = x1 - x2;         % differences
    sdDeltas = nanstd(deltas);  % standard deviation of the diffferences
    s        = sdDeltas;        % re-name
    d        =  meanDiff / s;   % Cohen's d (paired version)
    
end