function X_lim = chi_sq_limit (X, alpha)
% CHI_SQ_LIMIT One-tail confidence limit of X based on a chi squared distribution
%	with matching moments
%
% Syntax:
%	X_lim = chi_sq_limit (X, alpha)
%
% Inputs:
%	X:			Array of observations of the statistic X
%	alpha:		Significance level of the confidence limit
%
% Outputs:
%	X_lim:		One-tail confidence limit of X at signifcance level alpha
%

%{
chi_sq_limit.m
Version: 1.0.0
Date: 2023-12-19
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

chi_sq_limit.m
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

% Number of observations
N = size(X, 1);

% Mean of the statistic
mu = sum(X, 1)/N;
% Variance of the statistic
sigma = sum(((X - mu).^2), 1)/(N - 1);

% Degrees of freedom of the distribution
DOF = 2*(mu.^2)./sigma;
% Scale factor of the distribution
scalef = mu./DOF;
% Confidence limit
X_lim = scalef.*chi2inv(1 - alpha, DOF);
end