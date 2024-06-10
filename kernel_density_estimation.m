function [pdf, x_evp, h] = kernel_density_estimation (x_obs, varargin)
% KERNEL_DENSITY_ESTIMATION Univariate kernel density estimation of probability
%	distribution function using the Gaussian kernel
%
% Syntax:
%	pdf = kernel_density_estimation (x_obs)
%	pdf = kernel_density_estimation (x_obs, 'Key', value)
%	[pdf, h] = kernel_density_estimation (_)
%
% Inputs:
%	x_obs:	Obervations to be used to compute the PDF
%
% Outputs:
%	pdf:	Probability density function
%	x_evp:	Evaluation points of the probability density function
%	h:		Estimated bandwidth of the Gaussin kernel
%
% Keys and Values:
%	'Evaluation_points': x_evp
%		Column vector of points at which to evaluate the PDF
%	'Bandwidth': ['scott' | 'silverman' | h]
%		Bandwidth of the Gaussian kernel: optimal values can be computed by
%		Scott's rule ('scott', defualt) or Silverman's rule ('silverman'), or a
%		custom value can be provided by the user (h, as numeric)
%
%
% NOTE
%	The PDF of the random variable X estimated by KDE given a vecotr x_obs of N
%	observations of the random variable is defined as:
%		p(x) = (1/(N*h))*sum_n[k((x - x_obs(n))/h)]
%	where k is the gaussian kernel function:
%		k(x) = exp(-0.5*((x - x_obs(n)).^2)/h^2)/sqrt(2*pi)
%	and h is the kernel bandwidth. The kerneòl bandwidth can be provided by the
%	user. Alternatively, it can estimated by Scott's rule:
%		h = 1.059*bws*N^(-0.2)
%	or by Silverman's rule:
%		h = 0.9*bws*N^(-0.2)
%	In both cases, bws is the bandwidth scale factor, computer as:
%		bws = min([std(x_obs), iqr(x_obs)])
%
% NOTE
%	By defult, the evaluation points are defined as a set of 301 euqally spaced
%	points extending 3 bandwidts below the minium value in x_min and 3 bandwidths
%	above the maximum value in x.

%{
kernel_density_estimation.m
Version: 1.0.0
Date: 2023-12-20
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

kernel_density_estimation.m
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
x_evp = [];
bw_est = 'scott';

% Development cycle
if ~isempty(varargin)
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
			case 'Evaluation_points'
				x_evp = varargin{i + 1};
				if sum(size(x_evp) ~= 1) > 1
					error(['The function is meant for univariate kernel ...'...
						'density estimation'])
				end
			case 'Bandwidth'
				bw_est = varargin{i + 1};
				if isnumeric(bw_est)
					h = bw_est;
					bw_est = 'custom';
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

%% Initial operations

% Ajust dimensions of vectors for efficiency
x_obs = x_obs(:);
x_evp = x_evp(:);

% Number of observations
N = size(x_obs, 1);

% Number of evaultation points
Ne = 301;
% Span of the evaluation region
cutoff = 3;

%% Bandwidth estimation

if ~strcmp(bw_est, 'custom')
	% Variance of observation points
	s = sum(((x_obs - sum(x_obs, 1)/N).^2), 1)/(N - 1);
	% Inter-quartile-range of the obsevation points
	s_alt = iqr(x_obs, 1)/1.349;
	% bandwidth scale
	bws = min([s, s_alt]);

	% Bandwidth estimation
	switch bw_est
		case 'scott'
			% Use Scott's rule
			h = 1.059*bws*N^(-0.2);
		case 'silverman'
			% Use Silverman's rule
			h = 0.9*bws*N^(-0.2);
	end
end

%% Kernel density estimation

% Generate evaluation points
if isempty(x_evp)
	x_evp = linspace(min(x_obs) - cutoff*h, max(x_obs) + cutoff*h, Ne)';
end

% Probability density function by KDE
pdf = sum(gaussian_kernel(x_evp, x_obs', h), 2)/(N*h);

end