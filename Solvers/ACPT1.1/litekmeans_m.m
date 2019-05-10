function label = litekmeans_m(X, k)
% Perform k-means clustering.
%   X: d x n data matrix
%   k: number of seeds
% Written by Michael Chen (sth4nth@gmail.com).
% Copyright (c) 2009, Michael Chen
% All rights reserved.

% Modified by Wenzhao Zhao

i=0;
n = size(X,2);

if k==1
    label=ones(1,n);
	if verLessThan('matlab', '8.1') == 0 % transposes depending on MATLAB version
		label = label';
	end
else

	last = 0;
    label = ceil(k*rand(1,n));  % random initialization
	while any(label ~= last)
		[~,~,label] = unique(label);   % remove empty clusters
		E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
		center = X*(E*spdiags(1./sum(E,1)',0,k,k));    % compute center of each cluster
		last = label;
		[~,label] = max(bsxfun(@minus,center'*X,0.5*sum(center.^2,1)')); % assign samples to the nearest centers
		if verLessThan('matlab', '8.1') == 0 % transposes depending on MATLAB version
			label = label';
		end
		i=i+1;
 
	   if i>200
		   break;
	   end
	end

end

