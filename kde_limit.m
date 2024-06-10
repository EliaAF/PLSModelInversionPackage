function X_lim = kde_limit (X, alpha, varargin)
% KDE_SQ_LIMIT Confidence limit of X based on kernel density estimation wih
%	Gaussian kernel
%
% Syntax:
%	X_lim = chi_sq_limit (X, alpha)
%	X_lim = chi_sq_limit (X, alpha, 'Key', value)
%
% Inputs:
%	X:			Array of observations of the statistic X
%	alpha:		Significance level of the confidence limit
%
% Outputs:
%	X_lim:		Confidence limit of X at signifcance level alpha
%
% Keys and Values:
%	'Test_type': ['one_tail | 'two_tail]
%		Wheter the confidence limit should be evaulated as a one-tail test
%		(defualt) or a two-tail test
%	'Evaluation_points': x_evp
%		Column vector of points at which to evaluate the density function from
%		kernel density estimation
%	'Bandwidth': ['scott' | 'silverman' | h]
%		Bandwidth of the Gaussian kernel: optimal values can be computed by
%		Scott's rule ('scott', defualt) or Silverman's rule ('silverman'), or a
%		custom value can be provided by the user (h, as numeric)
%
%
% NOTE
%	The confidence limit returned by the function can be evaulted based on a
%	one-tail test or on a two-tail test. In the former case, only the upper
%	confidence limit is returned (the lower limit is assumed to be 0), while in
%	the latter case a column vecotr containing bothlower and upper confidence
%	limits is returned.
%
% NOTE
%	The mfunction yield the exact values of the diagnostic statistic
%	correpsonding to the declared significance level. As the probability density
%	functions is obstained in its sample form, a simple look-up by linear
%	interpolation is used to find the exact X_lim.

%{
kde_limit.m
Version: 1.0.0
Date: 2023-12-20
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

kde_limit.m
Copyright (c) 2020- Elia Arnese Feffin

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see https://www.gnu.org/licenses/gpl-3.0.html.
%}

%% Optional arguments development

% Optionals initialization
test_cl = 'one_tail';
x_evp = [];
bw_est = 'scott';

% Development cycle
if ~isempty(varargin)
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
			case 'Test_type'
				test_cl = varargin{i + 1};
				if strcmp(test_cl, 'one_tail') || strcmp(test_cl, 'two_tail')
				else
					error('Supported test are one-tail and two-tail')
				end
			case 'Evaluation_points'
				x_evp = varargin{i + 1};
				if sum(size(x_evp) ~= 1) > 1
					error(['The function is meant for univariate kernel ...'...
						'density estimation'])
				end
			case 'Bandwidth'
				bw_est = varargin{i + 1};
				if isnumeric(bw_est)
				else
					if strcmp(bw_est, 'scott') || strcmp(bw_est, 'silverman')
					else
						error(['Supported methods for bandwidth estimation '...
							"are Scott's rule and Silverman's rule, or a "...
						'custom bandwidth can be provided'])
					end
				end
			otherwise
				error(['Key ' key ' undefined'])
		end
	end
end

%% Probability density function by KDE

% Check if evaluation points have been provided
if isempty(x_evp)
	% Probability density function by KDE
	[pdf, x_evp, ~] = kernel_density_estimation(X,...
		'Bandwidth', bw_est...
	);
else
	% Probability density function by KDE
	[pdf, x_evp, ~] = kernel_density_estimation(X,...
		'Evaluation_points', x_evp,...
		'Bandwidth', bw_est...
	);
end

% Cumulative distribution function by KDE (integral by trapezoidal rule)
cdf = [0; cumsum(0.5*(pdf(1:(end - 1)) + pdf(2:end)).*diff(x_evp), 1)];

%% Confidence limit on the CDF

% Decide which confidence limit to estiamte
switch test_cl
	case 'one_tail'
		% Limit from significance level
		lim = 1 - alpha;
		% Find fist CDF entry grater that limit
		idx = find(cdf > lim, 1, 'first');
		% Find the extact confidence limit
		X_lim = x_evp(idx - 1) + (x_evp(idx) - x_evp(idx - 1))*(lim - cdf(idx - 1))/(cdf(idx) - cdf(idx - 1));
	case 'two_tail'
		% Limits from significance level
		lim = [alpha/2; 1 - alpha/2];
		% Pre-allocate confidence limit
		X_lim = zeros(2, 1);
		% Find fist CDF entry grater that lower limit
		idx = find(cdf > lim(1), 1, 'first');
		% Find the extact lower confidence limit
		X_lim(1) = x_evp(idx - 1) + (x_evp(idx) - x_evp(idx - 1))*(lim(1) - cdf(idx - 1))/(cdf(idx) - cdf(idx - 1));
		% Find fist CDF entry grater that upper limit
		idx = find(cdf > lim(2), 1, 'first');
		% Find the extact upper confidence limit
		X_lim(2) = x_evp(idx - 1) + (x_evp(idx) - x_evp(idx - 1))*(lim(2) - cdf(idx - 1))/(cdf(idx) - cdf(idx - 1));
end

end