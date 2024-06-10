function [c, ceq, nabla_c, nabla_ceq] = pls_IBO_nonlinear_constraint_wrapper (z, T_sq_bnd, SRE_bnd, NS_nl_bnd, Phi, T_sq_lim, Psi, SRE_lim, y_des, OO, OL, Q_tilde, Q_tilde_inv_L, Lambda_inv, RMSE, t_val, N, y_free)
% PLS_IBO_NONLINEAR_CONSTRAINT_WRAPPER Wrapper function for the nonlinear hard
%	constraints for PLS model inversion by optimization
%
% Syntax:
%	[c, ceq] = pls_IBO_nonlinear_constraint_wrapper (z, T_sq_bnd, SRE_bnd, NS_nl_bnd, Phi, T_sq_lim, Psi, SRE_lim, y_des, OO, OL, Q_tilde, Q_tilde_inv_L, Lambda_inv, RMSE, t_val, N, y_free)
%	[c, cew, nabla_c, nabla_ceq] = pls_IBO_nonlinear_constraint_wrapper (z, T_sq_bnd, SRE_bnd, NS_nl_bnd, Phi, T_sq_lim, Psi, SRE_lim, y_des, OO, OL, Q_tilde, Q_tilde_inv_L, Lambda_inv, RMSE, t_val, N, y_free)
%
% Inputs:
%	t_des:			Desired score (optimization variable) as column vector
%	T_sq_bnd:		Logical indicator for constraint of confidence region of T^2
%					(conisder the constraint if true, neglect otherwhise)
%	SRE_bnd:		Logical indicator for constraint of confidence region of SRE
%					(conisder the constraint if true, neglect otherwhise)
%	NS_nl_bnd:		Logical indicator for constraint of confidence region of the
%					null space (conisder the constraint if true, neglect
%					otherwhise)
%	Phi:			Scale matrix to compute the T^2 statistic
%	T_sq_lim:		Confidence limit of T^2
%	Psi:			Scale matrix to compute the SRE statistic
%	SRE_lim:		Confidence limit of SRE
%	y_des:			Target output for model inversion as column vector
%	OO:				Optimization-to-output variables matrix
%	OL:				Optimization-to-latent variables matrix
%	Q_tilde:		Modified output loadings matrix of the PLS model
%	Q_tilde_inv_L:	Left generalized inverse of the modified output loadings
%					matrix of the PLS model
%	Lambda_inv:		Scale matrix to compute the leverage of the solution
%	RMSE:			Root-mean squared errors of the output variables in
%					calibration (corrected for degrees of freedom of the model)
%	t_val:			Value of a t-distributed variables at the desides
%					significance level according to the two-tail test
%	N:				Number of observations in the calibration dataset of the PLS
%					model
%	y_free:			Logical vector marking free output variables (no target set)
%
% Outputs:
%	c:				Values of the inequality constraint functions at z
%	ceq:			Values of the equality constraint functions at z
%	nabla_c:		Gradients of of the inequality constraint functions at z
%	nabla_ceq:		Gradients of of the equality constraint functions at z
%
%
% NOTE
%	This function can be used to set hard constraints on PLS model inversion by
%	optimization using either the scores (z = t) or the inputs (z = x) as
%	optimization variables. The function is called internally by the primary
%	function pls_inversion_IBO, and calls itslef functions
%	pls_IBO_nonlinear_constraint_SD and pls_IBO_nonlinear_constraint_NS. See the
%	help to those three functions for details on the input arguments.
%
% NOTE
%	In order to handle cases where some of the output variables do not have an
%	equality constraint, therefore they are ``free'' (saved for possible
%	inequality constraints), such outpus are neglected both the in formulation of
%	the null space and in the estimation of its uncertainty. If the output i is
%	``free'', set the i-th component of the y_free vector to true. The y_free
%	vector should be a column vector with V_Y components, where V_Y is the number
%	of output variables.
%
% NOTE
%	There is fundamentally no useful nonlinear equality constraint that can be
%	set in PLS model inversion. However, ceq and nable_ceq must be provided
%	anyway due to the construction of the fmincon function in MatLab. In this
%	function, they are always returned as empty arrays.
%
% NOTE
%	The function returns c as a row vector: if all contraints are specified, the
%	first and second component refer to the hard constraint on T^2 and SRE,
%	respectively, of the optimization solution, the following V_Y components
%	refer to the lower limit of the null space, and the last V_Y components
%	bound the upper limit (V_Y is the number of output variables of the PLS
%	model). The gradients of the constraints functions, nabla_c, are returned as
%	a columns in a matrix nabla_c, where the i-th columns is the gradient of the
%	i-th constraint. Components of c and columns of nabla_c are removed for
%	unspecified constraints.

%{
pls_IBO_nonlinear_constraint_wrapper.m
Version: 1.0.0
Date: 2024-01-23
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_IBO_nonlinear_constraint_wrapper.m
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

% Dimensionality of the optimization variable
Dz = length(z);
% Numbers of output variables
V_Y = length(y_des);

% Initialize array for constraint values and gradients
c = NaN(1, 2 + 2*V_Y);
temp = NaN(2 + 2*V_Y, Dz);

% Assign constraints
if T_sq_bnd
	[c(1), temp(1, :)] = pls_IBO_nonlinear_constraint_SD(z, Phi, T_sq_lim);
end
if SRE_bnd
	[c(2), temp(2, :)] = pls_IBO_nonlinear_constraint_SD(z, Psi, SRE_lim);
end
if NS_nl_bnd
	% Assign constraints and gradients
	[c(3:end), temp(3:end, :)] = pls_IBO_nonlinear_constraint_NS(z, y_des, OO, OL, Q_tilde, Q_tilde_inv_L, Lambda_inv, RMSE, t_val, N);
	aux = [y_free; y_free];
else
	aux = true(2*V_Y, 1);
end
% Remove uneccesary constraints and gradients
idx = [~T_sq_bnd; ~SRE_bnd; aux];
c(idx) = [];
temp(idx, :) = [];

% No nonlinear equality constraints
ceq = [];

% Check if gradients are required
if nargout > 2
	% Assign gradients
	nabla_c = temp';
	% No nonlinear equality constraints
	nabla_ceq = [];
end
end