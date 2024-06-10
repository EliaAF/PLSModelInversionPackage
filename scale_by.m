function X_s = scale_by (X, mu, sigma)
% SCALE_BY Data arrays scaling with specified means and standard deviations
%
% Syntax:
%	X_s = scale_by (X, mu, sigma)
%
% Inputs:
%	X:		Data array to be scaled
%	mu:		Means to be used for scaling (row vector)
%	sigma:	Standard deviations to be used for scaling (row vector)
%
% Outputs:
%	X_s:	Scaled data array

%{
scale_by.m
Version: 1.0.0
Date: 2020-12-02
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

scale_by.m
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

% Replace null standard deviations
sigma(sigma == 0) = inf;

% Mean centring and standard deviation scaling
X_s = (X - mu)./repmat(sigma, size(X, 1), 1);
end