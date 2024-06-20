function[y,varargout] = winsor(x,p)
% WINSOR     Winsorize a vector 
% INPUTS   : x - n*1 data vector
%            p - 2*1 vector of cut-off percentiles (left, right) 
% OUTPUTS  : y - winsorized x, n*1 vector
%            i - (optional) n*1 value-replaced-indicator vector
% NOTES    : Let p1 = prctile(x,p(1)), p2 = prctile(x,p(2)). (Note
%            that PRCTILE ignores NaN values). Then 
%            if x(i) < p1, y(i) = min(x(j) | x(j) >= p1)
%            if x(i) > p2, y(i) = max(x(j) | x(j) <= p2)
% EXAMPLE  : x = rand(10,1), y = winsor(x,[10 90])
% AUTHOR   : Dimitri Shvorob, dimitri.shvorob@vanderbilt.edu, 4/15/07
% Copyright (c) 2016, Dimitri Shvorob
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
%POSSIBILITY OF SUCH DAMAGE.


if ~isvector(x)
   error('Input argument "x" must be a vector')
end  
if nargin < 2
   error('Input argument "p" is undefined')
end 
if ~isvector(p)
   error('Input argument "p" must be a vector')
end  
if length(p) ~= 2
   error('Input argument "p" must be a 2*1 vector')
end  
if p(1) < 0 || p(1) > 100
   error('Left cut-off percentile is out of [0,100] range')
end  
if p(2) < 0 || p(2) > 100
   error('Right cut-off percentile is out of [0,100] range')
end  
if p(1) > p(2)
   error('Left cut-off percentile exceeds right cut-off percentile')
end  
p = prctile(x,p);
i1 = x < p(1); v1 = min(x(~i1));
i2 = x > p(2); v2 = max(x(~i2));
y = x;
y(i1) = v1;
y(i2) = v2;
if nargout > 1
   varargout(1) = {i1 | i2};
end   