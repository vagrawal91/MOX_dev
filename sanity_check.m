function k = sanity_check(Dx,cvarX,q,r)
% VALIDATE_COMPONENTS - Sanity check for component selection
%
% Inputs:
%   Dx    - Optimal components from SVD analysis
%   cvarX - Cumulative variance vector from SVD of X
%   q     - Number of response variables
%   r     - Maximum allowable components, min(p,q)
%
% Output:
%   k     - Number of components to use

% Adjust Dx if it exceeds number of responses
if Dx > q
    Dx = q;
    warning('sanity_check: Dx adjusted with Dx=q. Variance captured: %.2f', cvarX(Dx)*100)
    if ceil(cvarX(Dx)*100) < 95
        warning(['Size of q is not satisfactory to capture the sufficient variance in X.' ...
         ' Hence, P and therefore, DQ may not be reliable.'])
    end
end

% Components to use
k = Dx;

% Check fundamental constraint
if k>r
    error('sanity_check: Violated k <= min(p,q): k=%d > r=%d',k,r)
end
end
