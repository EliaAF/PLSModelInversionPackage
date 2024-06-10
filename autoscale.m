function [X_s, mu, sigma] = autoscale (X)
% AUTOSCALE Data arrays autoscaling
%	Mean centering and strandrd deviation scaling
%
% Syntax:
%	[X_s, mu, sigma] = autoscale (X)
%
% Inputs:
%	X:		Data array to be autoscaled
%
% Outputs:
%	X_s:	Autoscaled data array
%	mu:		Means of the variables in X (row vector)
%	sigma:	Standard deviations of the variables in X (row vector)

%{
autoscale.m
Version: 1.0.0
Date: 2022-01-19
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

autoscale.m
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

% Means
mu = mean(X);
% Standard deviations for output
sigma = std(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Find logical columns
% logical_cols = all(X == 0 | X == 1, 1);
% % Scaling factor of logical columns are sqrt of probability of 1s in columns
% sigma_lcols = sqrt(sum(X, 1)/size(X, 1));
% % Join standard deviations of numerical and logical variables
% sigma = sigma.*~logical_cols + sigma_lcols.*logical_cols;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard deviations for autoscaling
sigma_s = sigma;

% Find very small standard deviations
nrm = sqrt(trace(X'*X));
nrm(nrm == 0) = 1/eps;
% Indexes for replacement
repl_std = (sigma./nrm) < eps*3;

% Replace useless standard deviations
sigma(repl_std) = 0;
sigma_s(repl_std) = inf;

% Mean centring and standard deviation scaling
X_s = (X - mu)./repmat(sigma_s, size(X, 1), 1);
end