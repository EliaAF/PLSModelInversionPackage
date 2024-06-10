function [c, nabla_c] = pls_IBO_nonlinear_constraint_SD (z, Zeta, SD_lim)
% PLS_IBO_NONLINEAR_CONSTRAINT_SD Nonlinear hard constraint on the confidence
%	limit of statistics defined as squared distances for PLS model inversion by
%	optimization
%
% Syntax:
%	c = pls_IBO_nonlinear_constraint_SD (z, Psi, SD_lim)
%	[c, nabla_c] = pls_IBO_nonlinear_constraint_SD (z, Psi, SD_lim)
%
% Inputs:
%	z:			Optimization variable as column vector
%	Zeta:		Scale matrix to compute the squared distance statistic
%	SD_lim:		Confidence limit of the squared distance (constraint value)
%
% Outputs:
%	c:			Value of the constraint function at z
%	nabla_c:	Gradient of of the constraint function at z
%
%
% NOTE
%	This function can be used to set hard constraints on PLS model inversion by
%	optimization using either the scores (t_des) or the inputs (x_des) as
%	optimization variables. Note that all vectors are assumed to be column
%	vectors in the following discussion.
% CASE 1: OPTIMIZATION IN THE SPACE OF LATENT VARIABLES
%	If the scores t are used as optimization variables, the solution is
%	restricted to the space of latent variables, therefore only T^2 can be
%	constrained. This can be done by setting:
%		z = t
%		Zeta = Sigma_inv = diag(sigma_sq.^-1)
%		SD_lim = T_sq_lim
%	where sigma_sq represent the vector of variances of the latent variables in
%	the calibration dataset, while T_sq_lim is the confidence limit of T^2. Note
%	that the squared reconstriction error (SRE) is forced to be null by
%	construction in this case.
% CASE 2: OPTIMIZATION IN THE INPUT SPACE
%	If the inputs x are used as optimization variables, the solution is allowed
%	to vary in the whole input space, therefore both T^2 and SRE can be
%	constrained. In the former case, set:
%		z = x
%		Zeta = Phi = W_ast*diag(sigma_sq.^-1)*W_ast'
%		SD_lim = T_sq_lim
%	where W_ast is the matrix of corrected input weights; in the latter case:
%		z = x_des
%		Zeta = Psi = (eye(V_X) - P*W_ast')'*(eye(V_X) - P*W_ast')
%		SD_lim = SRE_lim
%	where V_X is the number of input variables, P is the matrix of input
%	loadings, and SRE_lim is the confidence limit of SRE.
%
% NOTE
%	The function returns c as a scalar. The gradient of the constraint function,
%	nambla_c, is returned as a row vector, by the convention followed here.

%{
pls_IBO_nonlinear_constraint_SD.m
Version: 1.0.0
Date: 2024-01-18
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_IBO_nonlinear_constraint_SD.m
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

% Value of the constraint function
temp = z'*Zeta;
c = temp*z - SD_lim;

% Check if gradient is required
if nargout > 1
	% Compute gradient
	nabla_c = 2*temp;
end
end