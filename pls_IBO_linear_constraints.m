function [Alpha, beta, Alpha_eq, beta_eq] = pls_IBO_linear_constraints (L, h_lb, h_ub, h_eq, varargin)
% PLS_IBO_LINEAR_CONSTRAINTS Linear hard constraints for PLS model inversion by
%	optimization
%
% Syntax:
%	[Alpha, beta, Alpha_eq, beta_eq] = pls_IBO_linear_constraints (L, h_lb, h_ub, h_eq)
%	[Alpha, beta, Alpha_eq, beta_eq] = pls_IBO_linear_constraints (L, h_lb, h_ub, h_eq, Chlc_lb, hlc_lb, Chlc_ub, hlc_ul, Chlc_eq, hlc_eq)
%
% Inputs:
%	L:			Transfer matrix of the constraints to be set
%	h_lb:		Lower bounds for estimated input|latent|output variables
%	h_ub:		Upper bounds for estimated input|latent|output variables
%	h_eq:		Equality constraints for estimated input|latent|output variables
%	Chlc_lb:	Matrix of linear combination coefficients for lower bounds on
%				linear combinations of estimated input|output variables
%	hlc_lb:		Lower bounds for linear combinations of estimated
%				input|latent|output variables
%	Chlc_ub:	Matrix of linear combination coefficients for upper bounds on
%				linear combinations of estimated input|latent|output
%	hlc_ub:		Upper bounds for linear combinations of estimated
%				input|latent|output variables
%	Chlc_eq:	Matrix of linear combination coefficients for equality
%				constraints on linear combinations of estimated
%				input|latent|output
%	hlc_eq:		Equality constraints for linear combinations of estimated
%				input|latent|output variables
%
% Outputs:
%	Alpha:		Matrix of inequality constraints
%	beta:		Vector of inequality constraints
%	Alpha_eq:	Matrix of equality constraints
%	beta_eq:	Vector of equality constraints
%
%
% NOTE
%	The vectors to be used as bounds and equality constraints must be
%	pre-processed in the same way as the input|output data array used for model
%	calibration (meaning, columns must be centred on the means of the calibration
%	dataset and scaled on the variances of the calibration dataset).
%
% NOTE
%	This function can be used to set constraints on PLS model inversion by
%	optimization (in either the score or input space) on the estimated input,
%	latent, and output variables therefore variables h can be x, t, or y. Matrix
%	L should be the appropriate matrix to ``transfer'' the optimization variables
%	to the space of variables to be constrained. If the optimization problem is
%	formulated based on variables z, the following holds true.
%		- To constrain the output variables (h = y), the transfer matrix is set
%		  to the optimization-to-output variables, L = OO, such that:
%			• OO = Q_tilde	if z = t
%			• OO = B'		if z = x
%		- To constrain the latent variables (h = t), the transfer matrix is set
%		  to the optimization-to-latent variables, L = OL, such that:
%			• OL = eye(A)	if z = t
%			• OL = W_ast'	if z = x
%		- To constrain the input variables (h = x), the transfer matrix is set
%		  to the optimization-to-input variables, L = OI, such that:
%			• OI = P		if z = t
%			• OI = eye(V_X)	if z = x
%	Note that all vectors are regarded as column vectors in the discussion that
%	follows.
%
% NOTE
%	This function computes matrices and vectors for linear inequality constraints
%	(bounds) and equality constraints on estimated input|latent|output variables
%	and their linear combinations for PLS model inversion by optimization:
%		t_des = argmin_t [J(z)]
%			S.T.
%			Alpha*z <= beta
%			Alpha_eq*z = beta_eq
%			lb <= z <= up
%			c(z) <= 0
%			c_eq(z) = 0
%	where z is the optimization variable. The values of the constraints must be
%	provided as column vectors with V_X|A|V_Y components, where V_X|A|V_Y is the
%	number of input|latent|output variables. If some variables are to be left
%	unconstrained, set corresponding component to NaN. For example, if
%	V_X|A|V_Y = 3 and variable H_1 and H_3 are to be lower bounded, while
%	variable H_2 does not have a lower bound, then:
%		h_lb = [h_1_lb; NaN; h_2_lb]
%	Concerning constraints on linear combinations of input|latent|output
%	variables and taking upper bounds as examples, rows of Chlc_ub contain the
%	linear combination coefficients to compute the relevant linear combination,
%	while corresponding rows of hlc_ub contain the values of the constraints.
%	For example, if the estimated input|latent|output is to be constrained to the
%	unit simplex in a three dimensional space, one has to set lower bounds for
%	the three variables:
%		h_lb = [0; 0; 0]
%	and an upper bound on their linear combination such that:
%		H_1 + H_2 + H_3 <= 1
%	which can be written in matrix form as:
%		[1, 1, 1]*[H_1; H_2; H_3] <= 1
%	therefore:
%		Chlc_ub = [1, 1, 1]
%		hlc_ub	= 1
%	All constraints are stacked togther in matrices Alpha and Alpha_eq (and in
%	the corresponding vectors). Concerning inequality constraints, the order from
%	top to bottom is as follows: lower bounds on input|latent|output variables;
%	upper bounds on input|latent|output variables; lower bounds on linear
%	combinations of input|latent|output variables; upper bounds on linear
%	combinations of input|latent|output variables. Regarding equality
%	constraints, input|latent|output vairables are first, followed by their
%	lienar combinations. For bounds of input|latent|output variables, only rows
%	corresponding to bounded variables are reported in the Alpha matrix and
%	corresponding vector. If no inequality constraint is specified, then Alpha
%	and beta are returned as empty arrays. Same goes for Alpha_eq and beta_eq
%	if no equality constraint is specified.
%
% NOTE
%	Particularly interesting contraints for the PLS model inversion problem are
%	the ones that either snap the solution to the null space (if it exists) or
%	bound it within its confidence region. Regarding the former case, one simply
%	has to impose equality constraints as y_eq = y_des. For the latter case, the
%	method by Facco et al. (2015) can be used to estimate the uncertainty of the
%	null space, t_des_cl, and propagate it to the space of output variables to
%	get y_cl. Then set y_lb = y_des - y_cl and y_ub = y_des + y_cl.
%
% NOTE
%	The functions does not perform any check on redundant or infeasible
%	constraints. It is up to the user to make sure that all constraints work
%	together and to avoid making the optimization problem unsolvable.

%{
pls_IBO_linear_constraints.m
Version: 1.0.0
Date: 2024-01-12
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_IBO_linear_constraints.m
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

% Check if lower bounds are provided
if isempty(h_lb)
	Alpha_lb = [];
	beta_lb = [];
else
	% Make sure that z_lb is a column vector
	h_lb = h_lb(:);
	% Determine which input|output variables have lower bounds
	idx_lb = ~isnan(h_lb);
	% Matrix for lower bounds on inputs|outputs
	Alpha_lb = - L(idx_lb, :);
	% Vector for lower bounds on inputs|outputs
	beta_lb = - h_lb(idx_lb);
end

% Check if lower bounds are provided
if isempty(h_ub)
	Alpha_ub = [];
	beta_ub = [];
else
	% Make sure that z_ub is a column vector
	h_ub = h_ub(:);
	% Determine which input|output variables have upper bounds
	idx_ub = ~isnan(h_ub);
	% Matrix for upper bounds on inputs|outputs
	Alpha_ub = L(idx_ub, :);
	% Vector for upper bounds on inputs|outputs
	beta_ub = h_ub(idx_ub);
end

% Concatenate inequality constraint matrices and vectors
Alpha = [Alpha_lb; Alpha_ub];
beta = [beta_lb; beta_ub];

% Check if equality constraints are provided
if isempty(h_eq)
	Alpha_eq = [];
	beta_eq = [];
else
	% Make sure that z_eq is a column vector
	h_eq = h_eq(:);
	% Determine which input|output variables have equality constraints
	idx_eq = ~isnan(h_eq);
	% Matrix for upper bounds on inputs|outputs
	Alpha_eq = L(idx_eq, :);
	% Vector for upper bounds on inputs|outputs
	beta_eq = h_eq(idx_eq);
end

% Check if constraints on linear combinations of outputs are required
if ~isempty(varargin)
	% Develop arguments (and make sure vectors are columns)
	Chlc_lb = varargin{1};
	hlc_lb = varargin{2}(:);
	Chlc_ub = varargin{3};
	hlc_ub = varargin{4}(:);
	Chlc_eq = varargin{5};
	hlc_eq = varargin{6}(:);

	% Check if lower bounds are provided
	if isempty(hlc_lb)
		Alphalc_lb = [];
		betalc_lb = [];
	else
		% Matrix for lower bounds on inputs|outputs
		Alphalc_lb = - Chlc_lb*L;
		% Vector for lower bounds on inputs|outputs
		betalc_lb = - hlc_lb;
	end
	
	% Check if lower bounds are provided
	if isempty(hlc_ub)
		Alphalc_ub = [];
		betalc_ub = [];
	else
		% Matrix for upper bounds on inputs|outputs
		Alphalc_ub = Chlc_ub*L;
		% Vector for upper bounds on inputs|outputs
		betalc_ub = hlc_ub;
	end
	
	% Extend inequality constraint matrices and vectors
	Alpha = [Alpha; Alphalc_lb; Alphalc_ub];
	beta = [beta; betalc_lb; betalc_ub];
	
	% Check if equality constraints are provided
	if isempty(hlc_eq)
		Alphalc_eq = [];
		betalc_eq = [];
	else
		% Matrix for upper bounds on inputs|outputs
		Alphalc_eq = Chlc_eq*L;
		% Vector for upper bounds on inputs|outputs
		betalc_eq = hlc_eq;
	end

	% Extend equality constraint matrices and vectors
	Alpha_eq = [Alpha_eq; Alphalc_eq];
	beta_eq = [beta_eq; betalc_eq];
end
end