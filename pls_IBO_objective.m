function [J, nabla_J] = pls_IBO_objective (z_des, H, f)
% PLS_IBO_OBJECTIVE Objective function for PLS model inversion by optimization
%
% Syntax:
%	J = pls_IBO_objective (z_des, H, f)
%	[J, nabla_J] = pls_IBO_objective (z_des, H, f)
%
% Inputs:
%	z:			Optimization variable as column vector
%	H:			Matrix of the quadratic objective function
%	f:			Vector of the quadratic objective function
%
% Outputs:
%	J:			Value of the objective function at t_des
%	nable_J:	Gradient of of the objective function at t_des
%
%
% NOTE
%	This function gives a compact representation of the quadratic objective
%	function to be used on PLS model inversion by optimization using either the
%	scores (t_des) or the inputs (x_des) as optimization variables. Note that all
%	vectors are assumed to be column vectors in the following discussion.
% CASE 1: OPTIMIZATION IN THE SPACE OF LATENT VARIABLES
%	If the scores t_des are used as optimization variables, the solution is
%	restricted to the space of latent variables, therefore the objective
%	function can be formulated including a soft constraint on the T^2 statistic
%	of the solution:
%		J = (y_des - Q_tilde*t_des)'*Gamma_1*(y_des - Q_tilde*t_des) +
%			+ gamma_2*(t_des'*diag(sigma_sq.^-1)*t_des)
%	where y_des is the inversion target (assumed to be scaled on means and
%	standard deviations of output variables in the calibration dataset), Q_tilde
%	is the matrix of modified output loadings, and sigma_sq represents the vector
%	of variances of the latent variables in the calibration dataset. Gamma_1 is
%	the weighting matrix for the objective on the output target, while gamma_2 is
%	the weight for the soft contraint of T^2. The function can be written as a
%	standard quadratic form:
%		J = 0.5*t_des'*H*t_des + f'*t_des + c
%	where:
%		H = 2*(Q_tilde'*Gamma_1*Q_tilde + gamma_2*diag(sigma_sq.^-1))
%		f = - 2*Q_tilde'*Gamma_1*y_des
%		c = y_des'*Gamma_1*y_des
%	and c is irrelevant for the optimization. The gradient of the function is:
%		nabla_J = t_des'*H + f'
%	and, by the convention followed here, is a row vector.
% CASE 2: OPTIMIZATION IN THE INPUT SPACE
%	If the inputs x_des are used as optimization variables, the solution is
%	allowed to vary in the whole input space, therefore the objective function
%	can be formulated including soft constraints on both the T^2 statistic and
%	the squared reconstruction error (SRE) of the solution:
%		J = (y_des - B'*x_des)'*Gamma_1*(y_des - B'*x_des) +
%			+ gamma_2*(x_des'*W_ast*diag(sigma_sq.^-1)*W_ast'*x_des) +
%			+ gamma_3*(x_des'*(eye(V_X) - P*W_ast')'*(eye(V_X) - P*W_ast')*x_des)
%	where y_des B is the matrix of (outer) PLS regression coefficients, W_ast is
%	the matrix of corrected input weights, P is the matrix of input loadings, and
%	V_X is the number of input variables. Gamma_1 and gamma_2 have the same
%	meaning as in the previous case, while gamma_3 is the weight for the soft
%	contraint of SRE. The function can be written as a standard quadratic form:
%		J = 0.5*x_des'*H*x_des + f'*x_des + c
%	where:
%		H = 2*(B*Gamma_1*B' + gamma_2*(W_ast*diag(sigma_sq.^-1)*W_ast') +
%			+ gamma_3*((eye(V_X) - P*W_ast')'*(eye(V_X) - P*W_ast'))
%		f = - 2*B*Gamma_1*y_des
%		c = y_des'*Gamma_1*y_des
%	and c is irrelevant for the optimization. The gradient of the function is:
%		nabla_J = x_des'*H + f'
%	and, by the convention followed here, is a row vector.

%{
pls_IBO_objective.m
Version: 1.0.0
Date: 2024-01-18
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_IBO_objective.m
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

% Value of the objective function
temp = z_des'*H;
J = (0.5*temp + f')*z_des;

% Check if gradient is required
if nargout > 1
	% Compute gradient
	nabla_J = temp + f';
end
end