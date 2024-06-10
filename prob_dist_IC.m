function P_out = prob_dist_IC (IC_in)
% PROB_DIST_IC Probability distribution of information criteria
%
% Syntax:
%	P = prob_dist_IC (IC)
%
% Inputs:
%	IC:	Vector of values of an information criterion
%
% Outputs:
%	P:	Vecotr of probability distribution of information criteria

%{
prob_dist_IC.m
Version: 1.0.0
Date: 2021-06-29
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

prob_dist_IC.m
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

%% Input assignments

IC = IC_in;

%% Probability distributions

% Minimum of IC
min_IC = min(IC, [], 1);
% Offsets wrt min of IC
Delta_IC = IC - min_IC;
% Relative likelihoods
L_r = exp( - Delta_IC/2);
% Normalise relative likelihoods
P = L_r./sum(L_r, 1);

%% Output assignments

P_out = P;

end