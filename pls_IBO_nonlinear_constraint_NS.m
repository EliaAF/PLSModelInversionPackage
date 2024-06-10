function [c, nabla_c] = pls_IBO_nonlinear_constraint_NS (z, y_des, OO, OL, Q_tilde, Q_tilde_inv_L, Lambda_inv, RMSE, t_val, N)
% PLS_IBO_NONLINEAR_CONSTRAINT_NS Nonlinear hard constraint on the confidence
%	interval of the null space (method by Palací-López et al., 2019) for PLS
%	model inversion by optimization
%
% Syntax:
%	c = pls_IBO_nonlinear_constraint_NS (z, y_des, OO, OL, Q_tilde, Q_tilde_inv_L, Lambda_inv, RMSE, t_val, N)
%	[c, nabla_c] = pls_IBO_nonlinear_constraint_NS (z, y_des, OO, OL, Q_tilde, Q_tilde_inv_L, Lambda_inv, RMSE, t_val, N)
%
% Inputs:
%	z:				Optimization variable as column vector
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
%
% Outputs:
%	c:				Values of the constraint functions at z
%	nabla_c:		Gradients of of the constraint functions at z
%
%
% NOTE
%	This function can be used to set hard constraints on PLS model inversion by
%	optimization using either the scores (z = t) or the inputs (z = x) as
%	optimization variables. Internally, the function maps the optimization
%	variables onto the latent variables space by appropriate definitions of the
%	optimization-to-output variables matrix and optimization-to-latent variables
%	matrix passed as inputs. If the scores t are used as optimization variables,
%	the solution is restricted to the space of latent variables, therefore the
%	constraint is applied to z = t directly. On the other hand, if the inputs
%	x are used as optimization variables, the solution is allowed to vary in the
%	whole input space, therefore the constraint is applied to the projection of
%	the solution, t = W_ast'*x = OL*z, onto the latent space of the PLS model.
%	Note that several ``low level'' arguments are rquired as input to minimize
%	the numebr of repetitive calculations performed internally by the function,
%	thus to increase computational efficiency of the solution of the optimization
%	problem. In particular:
%		- Q_tilde_inv_L = Q_tilde'*pinv(Q_tilde*Q_tilde'), where the MatLab
%		  command pinv is used to generalize the formulation to all possible
%		  relationships among the number of output variables and the number of
%		  latent variables
%		- Lambda_inv = diag(sigma_sq.^1)/(N - 1), where siagma_sq is the vector
%		  of variances of the latent variables in the calibration dataset
%		- RMSE = sqrt(SSE/(N - dof)), where SSE is a vector of sum of squared
%		  errors of the output variables in calibration and dof are the degrees
%		  of freedom of the PLS model
%		- t_val = tinv(1 - alpha/2, N - dof), where alpha is the significance
%		  value of the confidence limit of the null space
%
% NOTE
%	The function returns c as a vector: the first V_Y components refer to the
%	lower limit of the null space, and the last V_Y components bound the upper
%	limit (V_Y is the number of output variables of the PLS model). The gradients
%	of the constraints functions, nabla_c, are returned as a rows in a matrix
%	nabla_c, where the i-th row is the gradient of the i-th constraint.

%{
pls_IBO_nonlinear_constraint_NS.m
Version: 1.0.0
Date: 2024-01-23
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_IBO_nonlinear_constraint_NS.m
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

% Predicted score
t_des = OL*z;
% Reconstructed y_des
y_des_pred = Q_tilde*t_des;
% Residual model error
r_des = y_des - y_des_pred;
% Perturbated score
t_tilde = t_des + Q_tilde_inv_L*r_des;
% Leverages of predictions
temp = t_tilde'*Lambda_inv;
h_y_des_pred = temp*t_tilde;
% Function of the leverage
f_lev = sqrt(1 + 1/N + h_y_des_pred);
% Confidence limit for the prediction
y_des_pred_cl = t_val*RMSE.*f_lev;

% Pre-allocation
c = zeros(2*V_Y, 1);
yrec = Q_tilde*t_tilde;
% Value of the constraint function: lower confidence limit
c(1:V_Y) = - y_des_pred + (yrec - y_des_pred_cl);
% Value of the constraint function: upper confidence limit
c(V_Y + 1:end) = y_des_pred - (yrec + y_des_pred_cl);

% Check if gradient is required
if nargout > 1
	% Pre-allocation
	nabla_c = zeros(2*V_Y, Dz);
	% Compute compnents of the gradient
	short_term = OL - Q_tilde_inv_L*OO;
	mid_term = Q_tilde*short_term;
	long_term = t_val*RMSE*(temp/f_lev)*short_term;
	% Compute gradient: lower confidence limit
	nabla_c(1:V_Y, :) = - OO + (mid_term - long_term);
	% Compute gradient: upper confidence limit
	nabla_c(V_Y + 1:end, :) = OO - (mid_term + long_term);
end
end