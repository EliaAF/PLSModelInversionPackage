function Q = Q_statistic (E)
% Q_STATISTIC Q statistic or squared prediction error
%
% Syntax:
%	Q = Q_statistic (E)
%
% Inputs:
%	E:			Array of residuals to be used for Q computation
%
% Outputs:
%	Q:			Q statistics for observations (rows) in Q
%

%{
Q_statistic.m
Version: 1.0.0
Date: 2023-12-19
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

Q_statistic.m
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

% Q statistic
Q = sum(E.^2, 2);
end