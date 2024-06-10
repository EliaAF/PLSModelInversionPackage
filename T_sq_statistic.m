function [T_sq, sigma_sq] = T_sq_statistic (T, varargin)
% T_SQ_STATISTIC Hotelling's T squared statistic
%
% Syntax:
%	T_sq = T_sq_statistic (T)
%	T_sq = T_sq_statistic (T, sigma_sq)
%	[T_sq, sigma_sq] = T_sq_statistic (_)
%
% Inputs:
%	T:			Array of observations to be used for T squared computation
%	sigma_sq:	Scale factors for T squared computation
%
% Outputs:
%	T_sq:		T squared statistics for observations (rows) in T
%	sigma_sq:	Scale factors for T squared computation (variances of columns of T)
%
%
% NOTE
%	The scale factors sigma_sq should be provided as a row vector with as many
%	components as the columns (variables) in T. If no sigma_sq is provided, the
%	variances of the columns of T are used as such and returned by the function,
%	otherwhise the sigma_sq provided as input is returned.

%{
T_sq_statistic.m
Version: 1.0.0
Date: 2023-12-19
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

T_sq_statistic.m
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

% Development cycle
if ~isempty(varargin)
	% Scale factors for T squared
	sigma_sq = varargin{1};
else
	% Number of observations
	N = size(T, 1);
	% Scale factors for T squared
	sigma_sq = sum(((T - sum(T, 1)/N).^2), 1)/(N - 1);
end

% Hotelling's T squared statistic
T_sq = sum((T.^2)*diag(sigma_sq.^(- 1)), 2);
end