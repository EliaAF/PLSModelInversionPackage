function [x_des_out, t_des_out, fval_out, ef_out, opt_info_out, lambda_out] = pls_inversion_IBO (y_des_in, Q_tilde_in, P_in, W_ast_in, sigma_sq_in, gamma_1_in, gamma_2_in, gamma_3_in, varargin)
% PLS_INVERSION_IBO PLS model inversion by optimization
%
% Syntax:
%	x_des = pls_inversion_IBO (y_des, Q_tilde, P, W_ast, sigma_sq, gamma_1, gamma_2, gamma_3)
%	x_des = pls_inversion_IBO (y_des, Q_tilde, P, W_ast, sigma_sq, gamma_1, gamma_2, gamma_3, nonzero_SRE)
%	x_des = pls_inversion_IBO (y_des, Q_tilde, P, W_ast, sigma_sq, gamma_1, gamma_2, gamma_3, nonzero_SRE, LCS)
%	x_des = pls_inversion_IBO (y_des, Q_tilde, P, W_ast, sigma_sq, gamma_1, gamma_2, gamma_3, nonzero_SRE, LCS, NS_eq, NS_lin_bnd, alpha, SSE, dof, N)
%	x_des = pls_inversion_IBO (y_des, Q_tilde, P, W_ast, sigma_sq, gamma_1, gamma_2, gamma_3, nonzero_SRE, LCS, NS_eq, NS_lin_bnd, alpha, SSE, dof, N, opts_qp)
%	x_des = pls_inversion_IBO (y_des, Q_tilde, P, W_ast, sigma_sq, gamma_1, gamma_2, gamma_3, nonzero_SRE, LCS, NS_eq, NS_lin_bnd, alpha, SSE, dof, N, opts_qp, T_sq_bnd, T_sq_lim, SRE_bnd, SRE_lim, NS_nl_bnd)
%	x_des = pls_inversion_IBO (y_des, Q_tilde, P, W_ast, sigma_sq, gamma_1, gamma_2, gamma_3, nonzero_SRE, LCS, NS_eq, NS_lin_bnd, alpha, SSE, dof, N, opts_qp, T_sq_bnd, T_sq_lim, SRE_bnd, SRE_lim, NS_nl_bnd, opts_nlp)
%	[x_des, t_des] = pls_inversion_IBO (_)
%	[x_des, t_des, fval] = pls_inversion_IBO (_)
%	[x_des, t_des, fval, ef, opt_info] = pls_inversion_IBO (_)
%	[x_des, t_des, fval, ef, opt_info, lambda] = pls_inversion_IBO (_)
%
% Inputs:
%	y_des:			Response variables to be used as targets in inversion
%	Q_tilde:		Matrix of modified output loadings of the PLS model
%	P:				Matrix of input loadings of the PLS model
%	W_ast:			Matrix of corrected input weights of the PLS model
%	sigma_sq:		Vector of variances of the latent variables in the
%					calibration dataset
%	gamma_1:		Vector of weights for individual targets variables in y_des
%	gamma_2:		Weight (scalar) for the soft contraint on the T^2 of z_des
%	gamma_3:		Weight (scalar) for the soft contraint on the SRE of z_des
%	nonzero_SRE:	Logical flag to carry out the optimization in the space of
%					input variables rather than in the space of latent variables
%	LCS:			Linear constraint structure
%	NS_eq:			Equality constraint on the null space
%	NS_lin_bnd:		Linear inequality constraints (bounds) on the null space
%	alpha:			Significance level of the confidence region of the null space
%	SSE:			Sum of squared errors of the output variables in calibration
%	dof:			Degrees of freedom of the PLS model
%	N:				Number of observations in the calibration dataset of the PLS
%					model
%	opts_qp:		Option structure for the quadratic programming solver
%	T_sq_bnd:		Nonlinear inequality constraint (bound) on the T^2 of z_des
%	T_sq_lim:		Confidence limit of the T^2 statistic
%	SRE_bnd:		Nonlinear inequality constraint (bound) on the SRE of z_des
%	SRE_lim:		Confidence limit of the SRE statistic
%	NS_nl_bnd:		Nonlinear inequality constraints (bounds) on the null space
%	opts_nlp:		Option structure for the nonlinear programming solver
%
% Outputs:
%	x_des:			Input variable vector estimated by model inversion
%	t_des:			Score vector estimated by model inversion
%	fval:			Value of the objective function at t_des
%	ef:				Status of the numerical optimizer (if used) at t_des
%	opt_info:		Additional information on the numerical optimization
%	lambda:			Langrange multipliers at t_des
%
%
% NOTE
%	The data array to be used as output target must be pre-processed in the same
%	way as the output data array used for model calibration (meaning, columns
%	must be centred on the means of the calibration dataset and scaled on the
%	variances of the calibration dataset). The array of input variables is
%	returned scaled as well.
%
% NOTE
%	The inversion method implemented in this function is based on the simplified
%	PLS model commonly used in the literature:
%		T = X*W_ast
%		X_est = T*P'
%		Y_est = T*Q_tilde'
%	where matrix Q_tilde is obtained by the output loading matrix (Q) and inner
%	regression coefficients (b) of the general PLS model as:
%		Q_tilde = Q*diag(b)
%	The function expects matrix Q_tilde to be provided as input. Note that the
%	PLS outer regression coefficients are still formulated as:
%		B = W_ast*Q_tilde'
%
% NOTE
%	Note that y_des and x_des (and t_des) are expected and returned as row
%	vectors. However, the function internally works with column vectors,
%	coherently with the formalism of optimization frameworks. The next notes
%	describe the problem trating y_des and x_des, and t_des (and likewise
%	vectors) as column vectors. It is thus helpful to re-formulate the PLS model
%	mentioned in the previous note in random vector form and using column
%	vectors:
%		t = W_ast'*x
%		x_rec = P*t
%		y_est = Q_tilde*t
%	where x is a the vector of input variables with dimensionality V_X, t is the
%	vector of latent variables (scores) with dimensionality A, and y is the
%	vector of output variables with dimesionality V_Y. x_rec denoted the
%	reconstruction of vector x by the PLS input data model, and y_est is the
%	estimated vector of output variables by the PLS regression model. Coherently,
%	the T^2 statistic is defined as:
%		T^2(x) = t'*diag(sigma_sq.^-1)*t = t'*W_ast'*diag(sigma_sq.^-1)*W_ast*t
%	where sigma_sq is a vector of variances of the latente variables, and the
%	squared reconstruction error (SRE) is defined as:
%		SRE = x*(eye(V_X) - W_ast*P')*(eye(V_X) - P*W_ast')*x
%
% NOTE
%	This function implements a generalized formulation of PLS model inversion by
%	optimization. Given target output (y_des), the general optimization problem
%	can be formulated as:
%		z_des = argmin_z [J(z)]
%			S.T.
%			Alpha*z <= beta
%			Alpha_eq*z = beta_eq
%			lb <= z <= up
%			c(z) <= 0
%			c_eq(z) = 0
%	where z is the optimization variable. The optimization can be carried out in
%	the space of latent variables (z = t, nonzero_SRE = false), thus the
%	optimization problem identifies the scorss t_des and the inputs are estimated
%	by the PLS input data model as x_des = P*t_des. In this case, x_des comforms
%	to the correlation structure modelled by PLS, thuse SRE(x_des) = 0. Note that
%	the input arguments W_ast and g_3 are not used in this case, thus they can be
%	set to []. Alternatively, the optimization can be carried out in the input
%	space directly (z = x, nonzero_SRE = true), thus SRE(x_des) is not
%	necessarily null. In either case, the objective function is generally defined
%	as:
%		J(z) = (y_des - y_est(z))'*Gamma_1*(y_des - y_est(z)) +
%				+ gamma_2*z'*Phi*z +
%				+ gamma_3*z'*Psi*z
%	where the first term accounts for the output target, the second term is a
%	soft contraint of the T^2 statistic of the solution, and the third term is a
%	soft constraint on the SRE of the solution. Gamma_1 is a matrix, thw diagonal
%	of which contains gamma_1, a vector of size V_Y (the number of output
%	variables). The components of gamma_1 are the weights for the targets on
%	individual variables in y_des. gamma_2 and gamma_3 are is a scalars weighting
%	the soft contraints on T^2 and SRE, respectively. Some matrices need to be
%	defined before proceeding.
%		- The optimization-to-output variables, OO, such that:
%			• OO = Q_tilde	if z = t
%			• OO = B'		if z = x
%		- The optimization-to-latent variables, OL, such that:
%			• OL = eye(A)	if z = t
%			• OL = W_ast'	if z = x
%		- The optimization-to-input variables, OI, such that:
%			• OI = P		if z = t
%			• OI = eye(V_X)	if z = x
%		- The scale matrix for T^2:
%			Phi = OL'*diag(sigma_sq.^-1)*OL
%		- The scale matrix for SRE:
%			Psi = OI'*(eye(V_X) - W_ast*P')*(eye(V_X) - P*W_ast')*OI
% Note that J(t) is a quadratic form in z, therefore the objective function can
%	be written as a standard quadratic form:
%		J(z) = 0.5*z'*H*z + z'*z + c
%	with:
%		H = 2*(OO'*Gamma_1*OO + gamma_2*Phi + gamma_3*Psi)
%		f = - 2*OO'*Gamma_1*y_des
%		c = y_des'*Gamma_1*y_des
%	Being PLS a linear model and J a quadratic form, the problem is convex and
%	admits a global optimum. Note that not all output variables need to have a
%	target value. In this case, set the corresponding component in y_des to NaN.
%	The complexity of the optimization problem varys based on the constraints set
%	on the problem itself.
% CASE 1
%	If no constraint is stated, the inversion can be solved as an unconstrained
%	quadratic programming, which admits an analytical solution:
%		z_des = H\(-f)
%	Note that if z = t, gamma_1 = ones(V_Y, 1), and gamma_2 = 0, the solution
%	reads:
%		t_des = (Q_tilde'*Q_tilde)\*Q_tilde'*y_des
%	which is the solution to direct inversion of PLS models with A < V_Y. To
%	generalize the solution of IBO to cases where A >= V_Y, the MatLab command
%	pinv is used for the matrix inversion in the equations above.
% CASE 2
%	If only linear contraints are stated, the inversion can be solved as a
%	general quadratic programming, which requires a numerical solution. The
%	quadprog function of the Optimization Toolbox is used to solve the problem.
%	Linear constraints on (estimated) input, latent, and output variables (and
%	their linear combinations) must be passed throught the LCS input argument.
%	This is a structure with fields described in the function
%	pls_IBO_initialize_constraints and used by the function
%	pls_IBO_linear_constraints to yield Alpha, beta, Alpha_eq, and beta_eq. See
%	the help fields of the the functions mentioned above for details on how to
%	provide values for the constraints. The present function internally assigns
%	the correct matrices (OO, OL, or OI) to ensure that the constraints are
%	applied to the correct space. Note that, if none of the linear constraints
%	defined by the LCS stucture are to be used, the LCS structure can be passed
%	anyway as imput to this function, as long as all its fiedls and subfileds are
%	empty. Additinal linear constraints on the null space (if it exists) can be
%	specified with the input arguments NS_eq (t_des is snapped to the null space
%	if NS_eq = true) and NS_lin_bnd (t_des is bounded to vary in the null space
%	confidence region at significance level alpha estimated with the method by
%	Facco et al. (2015) if NS_lin_bnd = true). Note that, in case only some of
%	the output variables have a target, the ``free'' outputs are not considered
%	in the formulation of the null space and its confidence region. Options can
%	be passed to the quadprog function using the opts_qp argument, which must be
%	a structure produced with the optimoptions command. Set the input argument
%	opts_qp to [] to retain the default options.
% CASE 3
%	For the problem to be a general quadratic programming (CASE 2), no nonlinear
%	constraint must be specificed. However, one may with to bound the solution
%	within a nonlinear formulation of the null space confidence region, such as
%	the one proposed by Palcí-López et al. (2019), or within the validity region
%	of the PLS model, meaning to bound T^2(z) and/or SRE(z) to be lesser or equal
%	to their respective confidence limits. The former constraint can be stated
%	setting the input argument NS_nl_bdn = true, while the latter constraints can
%	be stated by setting T_sq_bnd = true and/or SRE_bnd = true. These three
%	contraints are incorporated into the optimization problem as nonlinear
%	inequality constraints, thus through the nonlinear function c(t) <= 0. This
%	turns the problem into a general nonlinear programming, which requires a
%	numerical solution. The fmincon function of the Optimization Toolbox is used
%	to solve the problem. The constraint function c is called as the
%	pls_IBO_nonlinear_constraint_wrapper function, which itselfs calls
%	the functions pls_IBO_nonlinear_constraint_NS and
%	pls_IBO_nonlinear_constraint_SD. See the help fields of these functions for
%	details. Options can be passed to the fmincon function using the opts_nlp
%	argument, which must be a structure produced with the optimoptions command.
%	Set the input argument opts_nlp to [] to retain the default options.
%
% NOTE
%	The null space uncertainty is estimated with reference to the direct
%	inversion solution, as this solution serves as anchor point for the subspace
%	of solutions (that is the null space). In the case some output variables do
%	not have a target value, the corresponding rows of Q_tilde are nulled. This
%	effectlively accounts only for the target variables with targets, as the
%	``free'' ones do not contribute to the computation of scores in direct
%	inversion. It is easy to show empirically that the ``reduced null space''
%	one gets with the modified Q_tilde matrix still crosses the DI solution of
%	an inversion with all outputs specified. However, note that this approach is
%	merely a tricky approximation, and that it is not recommended to account for
%	the null space (especially for its uncertainty) in the case some of the
%	output variables do not have a specified target.
%
% NOTE
%	The last four output arguments are returned even for CASE 1, but are actually
%	meaningful only for CASE 2 and CASE 3. They correspond to the output
%	arguments of the MatLab functions quadprog and fmincon. See the
%	documentations of such functions for details.

%{
pls_inversion_IBO.m
Version: 1.0.0
Date: 2024-01-23
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_inversion_IBO.m
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

y_des = y_des_in;
Q_tilde = Q_tilde_in;
P = P_in;
W_ast = W_ast_in;
sigma_sq = sigma_sq_in;
gamma_1 = gamma_1_in;
gamma_2 = gamma_2_in;
gamma_3 = gamma_3_in;

%% Optinal argument development

% Optionals initialization
nonzero_SRE = [];
LCS = [];
NS_eq = [];
NS_lin_bnd = [];
alpha = [];
SSE = [];
dof = [];
N = [];
T_sq_bnd = [];
T_sq_lim = [];
SRE_bnd = [];
SRE_lim = [];
NS_nl_bnd = [];

% Default options for quadprog
opts_qp = optimoptions('quadprog',...
	'Display', 'none',...
	'Algorithm', 'active-set'...
);

% Default options for fmincon
opts_nlp = optimoptions('fmincon',...
	'Display', 'none',...
	'Algorithm', 'sqp',...
	'SpecifyObjectiveGradient', true,...
	'SpecifyConstraintGradient', true,...
	'CheckGradients', false...
);

if nargin > 8
	% CASE 1: no constraint
	nonzero_SRE = varargin{1};
	if nargin > 9
		% CASE 2: linear constraints only
		LCS = varargin{2};
		if nargin > 10
			% CASE 2: linear constraints on the null space
			NS_eq = varargin{3};
			NS_lin_bnd = varargin{4};
			alpha = varargin{5};
			SSE = varargin{6};
			dof = varargin{7};
			N = varargin{8};
			if nargin > 16
				% CASE 2: options for quadprog
				if ~isempty(varargin{9})
					opts_qp = varargin{9};
				end
				if nargin > 17
					% CASE 3: nonlinear constraints
					T_sq_bnd = varargin{10};
					T_sq_lim = varargin{11};
					SRE_bnd = varargin{12};
					SRE_lim = varargin{13};
					NS_nl_bnd = varargin{14};
					if nargin > 22
						% CASE 3: options for fmincon
						if ~isempty(varargin{15})
							opts_nlp = varargin{15};
						end
					end
				end
			end
		end
	end
end

%% Preliminary operations

% Numbers of output variables and of latent variables
[V_Y, A] = size(Q_tilde);
% Number of input variables
V_X = size(P, 1);

% Make sure dimensions are correct
y_des = y_des(:);
gamma_1 = gamma_1(:);

% Find free outputs, meaning without a target value
y_free = isnan(y_des);
% Adjust target vector and weight vector
y_des(y_free) = 0;
gamma_1(y_free) = 0;

% Initialize empty logical arguments for constraints
if isempty(NS_eq)
	NS_eq = false;
end
if isempty(NS_lin_bnd)
	NS_lin_bnd = false;
end
if isempty(T_sq_bnd)
	T_sq_bnd = false;
end
if isempty(SRE_bnd)
	SRE_bnd = false;
end
if isempty(NS_nl_bnd)
	NS_nl_bnd = false;
end

% Check if a null space exists
if (A <= V_Y - sum(y_free)) && (NS_eq || NS_lin_bnd || NS_nl_bnd)
	warning([...
		'No constraint can be set on the null space if the null space does '...
		'not exist: any constraint set on the null space will be neglected'...
	])
	NS_eq = false;
	NS_lin_bnd = false;
	NS_nl_bnd = false;
end

%% Determine IBO domain and case

% Check the optimization domain
if ~isempty(nonzero_SRE)
	if nonzero_SRE
		% Optimization in the space of input variables
		IBO_domain = 'input';
	else
		% Optimization in the space of latent variables
		IBO_domain = 'latent';
	end
else
	% Optimization in the space of latent variables
	IBO_domain = 'latent';
end
% If optmizing on latent variables, make sure that ``unused'' inputs make sense
if strcmp(IBO_domain, 'latent')
	W_ast = zeros(size(P));
	gamma_3 = 0;
end

% Gather constraint indicators in a cell array for ease of writing
opt_args_qp = {...
	NS_eq
	NS_lin_bnd
};
opt_args_nlp = {...
	T_sq_bnd
	SRE_bnd
	NS_nl_bnd
};

% Check is LCS is provided but empty
if isempty(LCS)
	LCS_empty = true;
	LCS = pls_IBO_initialize_constraints();
else
	if all(structfun(@(b) isempty(b) || all(isnan(b)), LCS.outputs)) && all(structfun(@(b) isempty(b) || all(isnan(b)), LCS.inputs)) && all(structfun(@(b) isempty(b) || all(isnan(b)), LCS.scores))
		LCS_empty = true;
	else
		LCS_empty = false;
	end
end

% Check if any agrument for quadratic programmng is provided
if all(cellfun(@(b) isempty(b) || b == false, opt_args_qp))
	QP_empty = true;
else
	QP_empty = false;
end

% Check if any agrument for nonlinear programmng is provided
if all(cellfun(@(b) isempty(b) || b == false, opt_args_nlp))
	NLP_empty = true;
else
NLP_empty = false;
end

% Determine case for optimization algorithm
if LCS_empty && QP_empty && NLP_empty
	% Unconstrained quadratic programming: analytical solution
	IBO_case = 'unconstrained';
elseif (~LCS_empty || ~QP_empty) && NLP_empty
	% Quadratic programming with linear constraints: numerical solution
	IBO_case = 'quadratic';
else% ~NLP_empty
	% Nonlinear programming: numerical solution
	IBO_case = 'nonlinear';
end

% % Display IBO domain and case
% disp([IBO_case, ' on ', IBO_domain])

%% Matrices for generalized problem

% Weighting matrix for output variables
Gamma_1 = diag(gamma_1);

% Inverse of score variances
Sigma_inv = diag(sigma_sq.^-1);
% Reconstruction matrix of PLS model
Rho = eye(V_X) - P*W_ast';

% Assign matrices based on optimization domain
switch IBO_domain
	case 'latent'
		% Dimensionality of the optimization space
		V_z = A;
		% Optimization-to-output variables matrix
		OO = Q_tilde;
		% Optimization-to-latent variables matrix
		OL = eye(A);
		% Optimization-to-input variables matrix
		OI = P;
	case 'input'
		% Dimensionality of the optimization space
		V_z = V_X;
		% Optimization-to-output variables matrix
		OO = Q_tilde*W_ast';
		% Optimization-to-latent variables matrix
		OL = W_ast';
		% Optimization-to-input variables matrix
		OI = eye(V_X);
end

% Scale matrix for T^2
Phi = OL'*Sigma_inv*OL;
% Scale matrix for SRE
Psi = OI'*(Rho'*Rho)*OI;

% Matrices of the objective function as quadratic form
H = 2*(OO'*Gamma_1*OO + gamma_2*Phi + gamma_3*Psi);
f = - 2*OO'*Gamma_1*y_des;

% Fix symmetry of the hessian
H = (H + H')/2;

% Analytical solution of the unconstrained problem
z_0 = pinv(H)*(-f); % We regularize to avoid the case-based structure

%% Matrices for linear constraints

% Check is numerical solutions is required
if strcmp(IBO_case, 'quadratic') || strcmp(IBO_case, 'nonlinear')

	% Develop linear constraints on output
	if all(structfun(@(b) isempty(b), LCS.outputs))
		% No linear constraints on outputs
		Alpha_Y = [];
		beta_Y = [];
		Alpha_eq_Y = [];
		beta_eq_Y = [];
	else
		% Assemble matrices for linear constraints on outputs
		[Alpha_Y, beta_Y, Alpha_eq_Y, beta_eq_Y] = pls_IBO_linear_constraints(...
			OO,...
			LCS.outputs.lower_bounds',...
			LCS.outputs.upper_bounds',...
			LCS.outputs.equality',...
			LCS.outputs.Cc_lower_bounds, LCS.outputs.Cv_lower_bounds,...
			LCS.outputs.Cc_upper_bounds, LCS.outputs.Cv_upper_bounds,...
			LCS.outputs.Cc_equality, LCS.outputs.Cv_equality...
		);
	end

	% Develop linear constraints on inputs
	if all(structfun(@(b) isempty(b), LCS.inputs))
		% No linear constraints on inputs
		Alpha_X = [];
		beta_X = [];
		Alpha_eq_X = [];
		beta_eq_X = [];
	else
		% Assemble matrices for linear constraints on inputs
		[Alpha_X, beta_X, Alpha_eq_X, beta_eq_X] = pls_IBO_linear_constraints(...
			OI,...
			LCS.inputs.lower_bounds',...
			LCS.inputs.upper_bounds',...
			LCS.inputs.equality',...
			LCS.inputs.Cc_lower_bounds, LCS.inputs.Cv_lower_bounds',...
			LCS.inputs.Cc_upper_bounds, LCS.inputs.Cv_upper_bounds',...
			LCS.inputs.Cc_equality, LCS.inputs.Cv_equality'...
		);
	end

	% Develop linear constraints on scores
	if all(structfun(@(b) isempty(b), LCS.inputs))
		% No linear constraints on scores
		Alpha_T = [];
		beta_T = [];
		Alpha_eq_T = [];
		beta_eq_T = [];
	else
		% Assemble matrices for linear constraints on scores
		[Alpha_T, beta_T, Alpha_eq_T, beta_eq_T] = pls_IBO_linear_constraints(...
			OL,...
			LCS.scores.lower_bounds',...
			LCS.scores.upper_bounds',...
			LCS.scores.equality',...
			LCS.scores.Cc_lower_bounds, LCS.scores.Cv_lower_bounds',...
			LCS.scores.Cc_upper_bounds, LCS.scores.Cv_upper_bounds',...
			LCS.scores.Cc_equality, LCS.scores.Cv_equality'...
		);
	end

	% Check if null space constraints are required
	if NS_eq || NS_lin_bnd || NS_nl_bnd
		% Reduced optimization-to-output variables matrix
		OO_red = OO;
		OO_red(y_free, :) = 0;

		% Check if any inequality constriant on the null space is required
		if NS_lin_bnd || NS_nl_bnd
			% Reduced output loading matrix
			Q_tilde_red = Q_tilde;
			Q_tilde_red(y_free, :) = 0;
			% Left generalized inverse
			Q_tilde_red_inv_L = Q_tilde_red'*pinv(Q_tilde_red*Q_tilde_red');
			% Scale matrix for leverages
			Lambda_inv = Sigma_inv/(N - 1);
			% Mean-squared error of the calibration dataset (DOF-corrected)
			RMSE = sqrt(SSE(:)/(N - dof));
			% Value of the t variable
			t_val = tinv(1 - alpha/2, N - dof);
		% Check if nonlinear optimization is required
		elseif strcmp(IBO_case, 'nonlinear')
			% Assign meaningless arrays for unused nonlinear constraint function
			Q_tilde_red = NaN(V_Y, A);
			Q_tilde_red_inv_L = NaN(A, V_Y);
			Lambda_inv = NaN(A, A);
			RMSE = NaN(V_Y, 1);
			t_val = NaN;
			N = NaN;
		end

		% Check if equality constraints on the null space are required
		if NS_eq
			NS_t_eq = y_des;
			NS_t_eq(y_free) = NaN;
		else
			NS_t_eq = [];
		end
		% Check if linear inequality constraints on the null space are required
		if NS_lin_bnd
			% Compute DI (possibly reduced) solution
			t_des_DI = Q_tilde_red_inv_L*y_des;
			% Leverages of predictions
			h_y_des_DI = t_des_DI'*Lambda_inv*t_des_DI;
			% Confidence limit for the prediction
			y_des_DI_cl = t_val*RMSE.*sqrt(1 + 1/N + h_y_des_DI);

			% Bounds for constraints
			NS_t_lb = y_des - y_des_DI_cl;
			NS_t_ub = y_des + y_des_DI_cl;
			NS_t_lb(y_free) = NaN;
			NS_t_ub(y_free) = NaN;
		else
			NS_t_lb = [];
			NS_t_ub = [];
		end
		
		% Assemble matrices for linear constraints on null space
		[Alpha_NS, beta_NS, Alpha_eq_NS, beta_eq_NS] = pls_IBO_linear_constraints(...
			OO_red,...
			NS_t_lb,...
			NS_t_ub,...
			NS_t_eq...
		);
	else
		% No linear constraints on the null space
		Alpha_NS = [];
		beta_NS = [];
		Alpha_eq_NS = [];
		beta_eq_NS = [];
		% Check if nonlinear optimization is required
		if strcmp(IBO_case, 'nonlinear')
			% Assign meaningless arrays for unused nonlinear constraint function
			OO_red = NaN(V_Y, V_z);
			Q_tilde_red = NaN(V_Y, A);
			Q_tilde_red_inv_L = NaN(A, V_Y);
			Lambda_inv = NaN(A, A);
			RMSE = NaN(V_Y, 1);
			t_val = NaN;
			N = NaN;
		end
	end

	% Assemble matrices for linear constraints
	Alpha = [Alpha_Y; Alpha_X; Alpha_T; Alpha_NS];
	beta = [beta_Y; beta_X; beta_T; beta_NS];
	Alpha_eq = [Alpha_eq_Y; Alpha_eq_X; Alpha_eq_T; Alpha_eq_NS];
	beta_eq = [beta_eq_Y; beta_eq_X; beta_eq_T; beta_eq_NS];
	
	% Bounds on optimizzation variables incorporated as linear constraints
	lb = [];
	ub = [];
end

%% PLS model inversion

% Determine case for IBO on scores
switch IBO_case
	case 'unconstrained'
		% Unconstrained quadratic programming: analytical solution
		z_des = z_0;
		% Additional output argiments
		fval = 0;
		ef = 1;
		opt_info = [];
		lambda = [];
	case 'quadratic'
		% Quadratic programming with linear constraints: numerical solution
		[z_des, fval, ef, opt_info, lambda] = quadprog(H, f,...
			Alpha, beta,...
			Alpha_eq, beta_eq,...
			lb, ub,...
			z_0,...
			opts_qp...
		);
	case 'nonlinear'
		% Nonlinear programming: numerical solution
		[z_des, fval, ef, opt_info, lambda] = fmincon(...
			@(z) pls_IBO_objective(z, H, f),...
			z_0,...
			Alpha, beta,...
			Alpha_eq, beta_eq,...
			lb, ub,...
			@(z) pls_IBO_nonlinear_constraint_wrapper(z,...
				T_sq_bnd,...
				SRE_bnd,...
				NS_nl_bnd,...
				Phi,...
				T_sq_lim,...
				Psi,...
				SRE_lim,...
				y_des,...
				OO_red,...
				OL,...
				Q_tilde_red,...
				Q_tilde_red_inv_L,...
				Lambda_inv,...
				RMSE,...
				t_val,...
				N,...
				y_free...
			),...
			opts_nlp...
		);
end

% Assign vectors of designed scores and input variables
x_des = OI*z_des;
t_des = OL*z_des;

% Reshape vectors as rows
x_des = x_des';
t_des = t_des';

%% Output assignments

x_des_out = x_des;
if nargout > 1
	t_des_out = t_des;
end
if nargout > 2
	fval_out = fval;
end
if nargout > 3
	ef_out = ef;
	opt_info_out = opt_info;
end
if nargout > 5
	lambda_out = lambda;
end
end