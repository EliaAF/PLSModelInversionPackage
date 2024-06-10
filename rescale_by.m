function X = rescale_by (X_s, mu, sigma)
% RESCALE_BY Data arrays rescaling with specified means and standard deviations
%
% Syntax:
%	X = rescale_by (X_s, mu, sigma)
%
% Inputs:
%	X_s:	Scaled data array to be rescaled
%	mu:		Means to be used for rescaling (row vector)
%	sigma:	Standard deviations to be used for rescaling (row vector)
%
% Outputs:
%	X:		Rescaled data array

%{
rescale_by.m
Version: 1.0.0
Date: 2021-12-02
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

rescale_by.m
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

% Mean recentring and standard deviation rescaling
X = X_s.*repmat(sigma, size(X_s, 1), 1) + mu;
end