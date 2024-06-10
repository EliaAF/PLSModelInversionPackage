function T_sq_lim = F_limit (T, alpha, dof)
% F_LIMIT One-tail confidence limit of T squared based on a F distribution
%
% Syntax:
%	T_sq_lim = F_limit (T, alpha, dof)
%
% Inputs:
%	T:			Array of observations to be used for confidence limit computation
%	alpha:		Significance level of the confidence limit
%	dof:		Degrees of freedom of the model
%
% Outputs:
%	T_sq_lim:	One-tail confidence limit of T squared at signifcance level alpha
%
%
% NOTE
%	The confidence limit based on the F distribution should be computed only for
%	statistics similar to a Hotelling's T squared statistic

%{
F_limit.m
Version: 1.0.0
Date: 2023-12-19
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

F_limit.m
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

% Number of observations and numebr of latent variables
N = size(T, 1);

% Confidence limit of T squared
T_sq_lim = (dof*(N - 1)*(N + 1)/(N*(N - dof)))*finv(1 - alpha, dof, N - dof);
end