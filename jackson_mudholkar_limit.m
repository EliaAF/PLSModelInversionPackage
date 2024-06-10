function Q_lim = jackson_mudholkar_limit (E, alpha)
% JACKSON_MUDHOLKAR_LIMIT One-tail confidence limit of Q based on the
%	Jackson-Mudholkar distribution
%
% Syntax:
%	Q_lim = F_limit (E, alpha)
%
% Inputs:
%	E:			Array of residuals to be used for confidence limit computation
%	alpha:		Significance level of the confidence limit
%
% Outputs:
%	Q_lim:		One-tail confidence limit of Q at signifcance level alpha
%
%
% NOTE
%	The confidence limit based on the Jackson-Mudholkar distribution should be
%	computed only for statistics similar to a Q statistic

%{
jackson_mudholkar_limit.m
Version: 1.0.0
Date: 2023-12-19
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

jackson_mudholkar_limit.m
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
N = size(T, 1);
% Variance of the residuals
sigma = sum(((E - sum(E, 1)/N).^2), 1)/(N - 1);

% Variational paramters of the distribution
theta = zeros(1, 3);
for j = 1:3
	theta(j) = sum(sigma.^j);
end
h0 = 1 - (2*theta(1)*theta(3))/(3*(theta(2)^2));
% Confidence limit of Q
Q_lim = theta(1)*(...
	1 +...
	norminv(1 - alpha)*(h0^2)*sqrt(2*theta(2))/theta(1) +...
	h0*(h0 - 1)*theta(2)/(theta(1)^2)...
)^(1/h0);
end