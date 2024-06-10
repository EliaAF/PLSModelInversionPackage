function model_out = build_pls (X_in, Y_in, A_in, varargin)
% BUILD_PLS PLS model calibration
%
% Syntax:
%	model = build_pls (X, Y, A)
%	model = build_pls (X, Y, A, 'Key', value)
%
% Inputs:
%	X:		Array of predictor variables to be used to calibrate the model
%	Y:		Array of response variables to be used to calibrate the model
%	A:		Number of latent variables
%
% Outputs:
%	model:	PLS model structure
%
% Keys and Values:
%	'Preprocessing': ['none' | 'mean_centre' | 'autoscale']
%		Preprocessing method to be applied to the input data structures, can be:
%		no scaling at all ('none'); mean-centring only ('mean_centre');
%		mean-centring and scaling to unit variance ('autoscale', default); note
%		that the algorithm implemented by this function requires mean centred
%		data, therefore if the option 'none' is used the user is in charge of
%		providing mean centred data, otherwise an error is issued
%	'Algorithm': ['simpls' | 'nipals']
%		Algorithm to be used for PLS calibration, can be either
%		statistically-inspired modification to partial least squares (SIMPLS) or
%		non-linear iteratrive partial least squared (NIPaLS) (default = 'nipals')
%	'Tol': [tol]
%		Tolerance for convergence as 2-norm on relative scores variation
%		(NIPaLS algorithm only, default = 1e-15)
%	'MaxIter': [max_iter]
%		Maximum number of iterations allowed (NIPaLS algorithm only,
%		default = 100)
%	'Normalise': ['standard' | 'loadings']
%		Normalisation scheme to be adopted, either standard normalisation
%		according to the algorithm of choice or PCA-coherent normalisation,
%		meaning both the X loadings and Y loadings are normalised (default =
%		'standard'); see functions pls_calibration_simpls and
%		pls_calibration_nipals for more information on normalisation schemes
%	'ErrBasedOn': ['unscaled' | 'scaled']
%		Scaling of predictions, reconstructions and errors, whether they
%		should be returned and computed from scaled or unscaled entities (default
%		= 'unscaled'); note that SSE, RMSE, bias and SE will be reported as
%		scaled even in the unscaled errors are requested as they are meant to
%		assess the performance of the model and the comparison of errors on
%		different variables is possible only if they are all scaled
%	'Contrib': ['simple' | 'absolute']
%		Kind of contributions to diagnostics, which can be simple indictors
%		either positive or negative according to the direction of deviation and
%		with two-tail confidence limits (approach of Miller, 1998), or purely
%		positive contributions (like "squadred") that sum up the diagnostic
%		values and with one-tail confidence limits (approach of Westerhuis, 2000)
%		(default = 'simple')
%	'ConLim': [lim]
%		Confidence limit for diagnostics statistics (default = 0.95)
%	'DOFMethod:['naive' | 'cv_based' | 'krylov']
%		Method for computing degrees of freedom of the model, can be either based
%		on the number of latent variables ('naive', dof = A + 1, default), based
%		on leave-one-out cross-validation (cv_based, dof estimated according
%		to Van Der Voet, 1999) or based on Krylov sequence ('krylov', dof
%		estimated according to Kramer, 2012)
%	'TsqLimMethod': ['chisq' | 'F' | 'KDE']
%		Method for computing confidence limits on T_sq_T and T_sq_U, can be
%		either the chi squared distribution method (default), on the F
%		distribution method or on kernel density estimation
%	'SqErrLimMethod': ['chisq' | 'jack_mud' | 'KDE']
%		Method for computing confidence limits on SPE_Y, SRE_X and SRE_Y, can be
%		either the chi squared distribution method (default), on the
%		Jackson-Mudholkar equation or on kernel density estimation
%	'ContribLimMethod': ['norm' | 't']
%		Method for computing the confidence limits of contributions to
%		diagnostics, can be based on a normal distribution or on a t distribution
%		(default = 'norm')
%	'EllipseForPlots: ['two' | 'full']
%		Method for computing the semiaxes of the confidence ellipse for scores
%		plot, compute from a F distribution with considering only two latent
%		variables (A_F = 2) or all the requested ones (A_F = A, default)
%	'ObsNames': [obs_names]
%		Names of the observations as chars in a cell array (default are
%		progressive numerical identifiers prefixed by the letter O)
%	'XVarNames': [X_var_names]
%		Names of the predictor variables as chars in a cell array (default are
%		progressive numerical identifiers prefixed by the letter X)
%	'YVarNames': [Y_var_names]
%		Names of the response variables as chars in a cell array (default are
%		progressive numerical identifiers prefixed by the letter Y)
%
%
% NOTE
%	PLS offers models for both X and Y:
%		X = T*P' + E
%		Y = U*Q' + F
%	as well as a model for their interrelation:
%		Y = X*B + R
%	where:
%		E is the residual matrix for the decomposition of X
%		F is the residual matrix for the decomposition of Y
%		R is the residual matrix for the regression of Y from X
%	The workflow for prediction is a loop-based procedure:
%		1) Initialise the predicition y_approx = 0
%		2) For each principal component a = 1:A:
%			a) Compute the score using the weight a:
%				t_a = x*w_a;
%			b) Deflate the sample using the loading a:
%				x = x - t_a*p_a'
%			c) Regress the Y scores from the X scores by means of the
%			   inner-relation regression coefficients and compose the
%			   contribution a to the predicition using the Y loadings,
%			   accumulating it on y_approx:
%				y_approx = y_approx + t_a*b_a*q_a';
%	The function also implement the modified nipals algorithm, which returns the
%	corrected weights. In such a case, the prediction workflow can be summarised
%	into a single operation defining the PLS regression coefficients matrix:
%		B = W_ast*diag(b)*Q'
%	thus:
%		y_approx = x*B
%
% NOTE
%	This function also serves as fundation for building PLSDA models. For such
%	occasions, the Y matrix should be prepared using the class_coding function; a
%	PLS model is then to be built using the present function and its output
%	passed to the pls_to_plsda function fo conversion into a PLSDA model.

%{
build_pls.m
Version: 1.0.1
Date: 2022-12-20
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

build_pls.m
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

X_unscaled = X_in;
Y_unscaled = Y_in;
A = A_in;

%% Initial checks

% Check if the data array is unfolded
if size(X_unscaled, 3) ~= 1 || size(Y_unscaled, 3) ~= 1
	error('The data arrays must be unfolded for PLS model calibration')
end
% Check if the requested number of LVs is feasible
if A > min(size(X_unscaled))
	error(['The number of latent variables cannot exceed min(size(X)) = '...
		num2str(min(size(X_unscaled)))])
end
% Check if there are missing values
if sum(ismissing(X_unscaled), 'all') ~= 0 || sum(ismissing(Y_unscaled), 'all') ~= 0
	error('Missing values found: PLS cannot be calibrated')
end

% Number of observations and number of variables
[N, V_X] = size(X_unscaled);
V_Y = size(Y_unscaled, 2);

%% Optional arguments development

% Optionals initialization
preprocess = 'autoscale';
alg = 'nipals';
tol = 1e-15;
max_iter = 100;
norm_sch = 'standard';
err_on = 'unscaled';
contrib = 'simple';
lim = 0.95;
dof_method = 'naive';
Tsq_lim_method = 'chisq';
SqE_lim_method = 'chisq';
con_lim_method = 'norm';
l_kind = 'full';

obs_names = cellstr([repmat('O', N, 1) char(pad(replace(string(num2str((1:N)')),  ' ', ''), length(num2str(N)), 'left', '0'))]);
X_var_names = cellstr([repmat('X', V_X, 1) char(pad(replace(string(num2str((1:V_X)')),  ' ', ''), length(num2str(V_X)), 'left', '0'))]);
Y_var_names = cellstr([repmat('Y', V_Y, 1) char(pad(replace(string(num2str((1:V_Y)')),  ' ', ''), length(num2str(V_Y)), 'left', '0'))]);

% Development cycle
if ~isempty(varargin)
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
			case 'Preprocessing'
				preprocess = varargin{i + 1};
				if strcmp(preprocess, 'none') || strcmp(preprocess, 'mean_centre') || strcmp(preprocess, 'autoscale')
				else
					error(['Supported preprocessing methods are no scaling, '...
						'mean-centring and autoscaling'])
				end
			case 'Algorithm'
				alg = varargin{i + 1};
				if strcmp(alg, 'simpls') || strcmp(alg, 'nipals')
				else
					error('Supported algorithms for PLS are SIMPLS and NIPaLS')
				end
			case 'Tol'
				tol = varargin{i + 1};
			case 'MaxIter'
				max_iter = varargin{i + 1};
			case 'Normalise'
				norm_sch = varargin{i + 1};
				if strcmp(norm_sch, 'standard') || strcmp(norm_sch, 'loadings')
				else
					error('Supported normalisations schemes are standard and loadings')
				end
			case 'ErrBasedOn'
				err_on = varargin{i + 1};
				if strcmp(err_on, 'scaled') || strcmp(err_on, 'unscaled')
				else
					error('Undefined key for ErrBasedOn')
				end
			case 'Contrib'
				contrib = varargin{i + 1};
				if strcmp(contrib, 'simple') || strcmp(contrib, 'absolute')
				else
					error(['Supported kinds of contributions to diagnostics '...
						'are simple and absolute'])
				end
			case 'ConLim'
				lim = varargin{i + 1};
			case 'DOFMethod'
				dof_method = varargin{i + 1};
				if strcmp(dof_method, 'naive') || strcmp(dof_method, 'cv_based') || strcmp(dof_method, 'krylov')
				else
					error(['Supported methods for degrees of freedom '...
						'estimation are naive, cv_based and krylov'])
				end
			case 'TsqLimMethod'
				Tsq_lim_method = varargin{i + 1};
				if strcmp(Tsq_lim_method, 'chisq') || strcmp(Tsq_lim_method, 'F') || strcmp(Tsq_lim_method, 'KDE')
				else
					error(['Supported methods for limits on T_sq_T and '...
						'T_sq_U are chi^2 method, F distribution method '...
						'and kernel density estimation method'])
				end
			case 'SqErrLimMethod'
				SqE_lim_method = varargin{i + 1};
				if strcmp(SqE_lim_method, 'chisq') || strcmp(SqE_lim_method, 'jack_mud') || strcmp(SqE_lim_method, 'KDE')
				else
					error(['Supported methods for limits on SPE_Y, SRE_X and '...
						'SRE_Y are chi^2 method, Jackson-Mudholkar equation '...
						'and kernel density estimation method'])
				end
			case 'ContribLimMethod'
				con_lim_method = varargin{i + 1};
				if strcmp(con_lim_method, 'norm') || strcmp(con_lim_method, 't')
				else
					error(['Supported methods for contribution limits are '...
						'normal distribution method or t distribution method'])
				end
			case 'EllipseForPlots'
				l_kind = varargin{i + 1};
				if strcmp(l_kind, 'two') || strcmp(l_kind, 'full')
				else
					error('Unrecognised value of the EllipseForPlots key')
				end
			case 'ObsNames'
				obs_names = varargin{i + 1};
				if any(size(obs_names) == N)
				else
					error(['Number of observation labels does not match the '...
						'number of observations'])
				end
			case 'XVarNames'
				X_var_names = varargin{i + 1};
				if any(size(X_var_names) == V_X)
				else
					error(['Number of predictor variable labels does not '...
						'match the number of X variables'])
				end
			case 'YVarNames'
				Y_var_names = varargin{i + 1};
				if any(size(Y_var_names) == V_Y)
				else
					error(['Number of response variable labels does not '...
						'match the number of Y variables'])
				end
			otherwise
				error(['Key ' key ' undefined'])
		end
	end
end

% Discourage on using cv_based DOF for multivariate Y
if strcmp(dof_method, 'cv_based') && V_Y > 1
	warning(['cv_based DOF estimation is meant for univariate Y: a variation '...
		'for multivariate Y will be used. However, for multivariate Y it is '...
		'suggested to use naive DOF estimation'])
end

%% Model structure initialization

% Model structure
model = initialize_pls;

% Dimesions and labels assignments
model.dimensions.N = N;
model.dimensions.V_X = V_X;
model.dimensions.V_Y = V_Y;
model.dimensions.A = A;
model.info.obs_names = obs_names;
model.info.X_var_names = X_var_names;
model.info.Y_var_names = Y_var_names;
model.info.preprocessing = preprocess;
model.info.algorithm = alg;
model.info.normalisation = norm_sch;
model.info.error_based_on = err_on;
model.info.contribution_method = contrib;
model.info.dof_method = dof_method;
model.info.Tsq_lim_method = Tsq_lim_method;
model.info.SqE_lim_method = SqE_lim_method;
model.info.con_lim_method = con_lim_method;
model.info.l_kind = l_kind;

%% Preprocessing

% Choiche of the preprocessing method
switch preprocess
	case 'none'
		mu_X = zeros(1, V_X);
		sigma_X = ones(1, V_X);
		X = X_unscaled;
		mu_Y = zeros(1, V_Y);
		sigma_Y = ones(1, V_Y);
		Y = Y_unscaled;
	case 'mean_centre'
		mu_X = mean(X_unscaled);
		sigma_X = ones(1, V_X);
		X = X_unscaled - mu_X;
		mu_Y = mean(Y_unscaled);
		sigma_Y = ones(1, V_Y);
		Y = Y_unscaled - mu_Y;
	case 'autoscale'
		[X, mu_X, sigma_X] = autoscale(X_unscaled);
		[Y, mu_Y, sigma_Y] = autoscale(Y_unscaled);
end

% Assigments to the model structure
model.data.X = X;
model.data.Y = Y;
model.data.X_uns = X_unscaled;
model.data.Y_uns = Y_unscaled;
model.scaling.mu_X = mu_X;
model.scaling.sigma_X = sigma_X;
model.scaling.mu_Y = mu_Y;
model.scaling.sigma_Y = sigma_Y;

%% Model calibration

% Selection of the algorithm for model calibration
switch alg
	case 'simpls'
		[B, W_ast, P, T, Q, U, b] = pls_calibration_simpls (X, Y, A, 'Normalise', norm_sch);
		W = [];
	case 'nipals'
		[B, W_ast, P, T, Q, U, b, W] = pls_calibration_nipals(X, Y, A, 'Tol', tol, 'MaxIter', max_iter, 'Normalise', norm_sch);
end

% Variances of X and Y scores for diagnostics
sigma_sq_T = sum(((T - sum(T, 1)/N).^2), 1)/(N - 1);
sigma_sq_U = sum(((U - sum(U, 1)/N).^2), 1)/(N - 1);

% Assigments to the model structure
model.parameters.B = B;
model.parameters.W_ast = W_ast;
model.parameters.P = P;
model.parameters.b = b;
model.parameters.Q = Q;
model.parameters.W = W;
model.parameters.sigma_sq_T = sigma_sq_T;
model.parameters.sigma_sq_U = sigma_sq_U;
model.prediction.T = T;
model.prediction.U = U;

%% Model performance

% Pre-allocation
PV_Y = zeros(A, 1);
EV_X = zeros(A, 1);
EV_Y = zeros(A, 1);
R_sq_Y = zeros(A, V_Y);

SSE = zeros(A, V_Y);
bias = zeros(A, V_Y);
SE = zeros(A, V_Y);

% Predicted Y
Y_pred = X*B;
% Reconstructued matrices
X_rec = T*P';
Y_rec = U*Q';
% Y prediction residuals
R = Y - Y_pred;
% Reconstruction errors
E = X - X_rec;
F = Y - Y_rec;

% Total variances in X and Y for EVs
ssq_X = trace(X'*X);
ssq_Y = trace(Y'*Y);
% Total variance for each Y for R_sq
ssq_Y_pv = diag(Y'*Y)';
% Loop over latent variables
for a = 1:A
	% Prediction error with LV a
	R_a = Y - T(:, a)*b(a)*Q(:, a)';
	% Reconstruction errors with LV a
	E_a = X - T(:, a)*P(:, a)';
	F_a = Y - U(:, a)*Q(:, a)';
	% Explained variance in prediction of Y
	PV_Y(a) = 1 - trace(R_a'*R_a)/ssq_Y;
	% Explained variances in reconstruction
	EV_X(a) = 1 - trace(E_a'*E_a)/ssq_X;
	EV_Y(a) = 1 - trace(F_a'*F_a)/ssq_Y;
	% Determination coefficient on each response variable
	R_sq_Y(a, :) = 1 - sum(R_a.^2, 1)./ssq_Y_pv;
	% Prediction error with LVs up to a
	R_a = Y - T(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
	% SSE in calibration (prediction error on Y)
	SSE(a, :) = sum(R_a.^2, 1);
	% Error bias
	bias(a, :) = sum(R_a, 1)/N;
	% Standard error
	SE(a, :) = sqrt(sum((R_a - bias(a, :)).^2, 1)/(N - 1));
end
% Cumulative explained variance in prediction of Y
CPV_Y = cumsum(PV_Y);
% Cumulative explained variances in reconstruction
CEV_X = cumsum(EV_X);
CEV_Y = cumsum(EV_Y);
% Cumulative determination coefficient on each response variable
CR_sq_Y = cumsum(R_sq_Y);
% Root mean squared error
RMSE = sqrt(SSE/N);
% Extracted eigenvalues by T and U (non-normalised EV_X and PV_Y)
switch alg
	case 'simpls'
		switch norm_sch
			case 'standard'
				lambda_T = diag(P'*P)/(N - 1);
				lambda_U = diag(Q'*Q)/(N - 1);
			case 'loadings'
				lambda_T = diag(T'*T)/(N - 1);
				lambda_U = diag(U'*U)/(N - 1);
		end
	case 'nipals'
		switch norm_sch
			case 'standard'
				lambda_T = ssq_X*EV_X/(N - 1);
				lambda_U = ssq_Y*EV_Y/(N - 1);
			case 'loadings'
				lambda_T = diag(T'*T)/(N - 1);
				lambda_U = diag(U'*U)/(N - 1);
		end
end
% Variance extracted by eigenvalues on Y appears to be linked to different paths
% in different algorithms. Look at the equations for NIPaLS, standard
% normalisation scheme: in NIPaLS eigenvalues of Y are related to the explained
% variance in reconstruction, while in SIMPLS that equation does not work: EV_Y
% have to be replaced with PV_Y, the explained variance in prediction. It appears
% that SIMPLS is actually aimed at prediction as U has no meaning at all in that
% algorithm. In fact it is also impossible to compute EV_Y with the general
% equation, unless U is normalised first.

% Copy entites for diagnostics
R_fd = R;
E_fd = E;
F_fd = F;

% Rescale entities according to preferences
if strcmp(err_on, 'unscaled')
	Y_pred = Y_pred*diag(sigma_Y) + mu_Y;
	X_rec = X_rec*diag(sigma_X) + mu_X;
	Y_rec = Y_rec*diag(sigma_Y) + mu_Y;
	R = R*diag(sigma_Y);
	E = E*diag(sigma_X);
	F = F*diag(sigma_Y);
end

% Assigments to the model structure
model.prediction.Y_pred = Y_pred;
model.prediction.X_rec = X_rec;
model.prediction.Y_rec = Y_rec;
model.prediction.R = R;
model.prediction.E = E;
model.prediction.F = F;
model.performance.PV_Y = PV_Y;
model.performance.CPV_Y = CPV_Y;
model.performance.EV_X = EV_X;
model.performance.CEV_X = CEV_X;
model.performance.EV_Y = EV_Y;
model.performance.CEV_Y = CEV_Y;
model.performance.R_sq_Y = R_sq_Y;
model.performance.CR_sq_Y = CR_sq_Y;
model.performance.SSE = SSE;
model.performance.RMSE = RMSE;
model.performance.bias = bias;
model.performance.SE = SE;
model.performance.lambda_T = lambda_T;
model.performance.lambda_U = lambda_U;

%% Model selection criteria

% Degrees of freedom of the model
[dof, ~] = pls_dof(X_unscaled, Y_unscaled, A,...
	'Method', dof_method,...
	'Preprocessing', preprocess,...
	'Algorithm', alg,...
	'Tol', tol,...
	'MaxIter', max_iter,...
	'Normalise', norm_sch...
);

% Adjusted determination coefficient on each response variable
R_sq_adj_Y = 1 - ((N - 1)./(N - dof)).*(1 - R_sq_Y);

% Information criteria
[IC, pIC] = information_criteria(SSE, dof, N);

% Assigments to the model structure
model.selection.dof = dof;
model.selection.R_sq_adj_Y = R_sq_adj_Y;
model.selection.IC = IC;
model.selection.pIC = pIC;

% Compress DOF for multivariate Y
dof = sum(dof(end, :), 2)/V_Y;

%% Model diagnostics

% Pre-allocation
VIP = zeros(V_X, V_Y);

% Hotelling's T^2 of T and U
T_sq_T = T_sq_statistic(T, sigma_sq_T);
T_sq_U = T_sq_statistic(U, sigma_sq_U);
% Squared prediction error (SPE) on Y
SPE_Y = Q_statistic(R_fd);
% Squared reconstruction errors (SRE) of X and Y
SRE_X = Q_statistic(E_fd);
SRE_Y = Q_statistic(F_fd);

% Choice of the kind of contributions
switch contrib
	case 'simple'
		% Upper and lower confidence limits
		lim_con = lim + (1 - lim)/2;
		% Contributions to Hotelling's T^2 for T and U: reconstruction path
		T_sq_T_con = (...
			T*sqrt(diag(sigma_sq_T.^(- 1)))...
		)*(P*pinv(P'*P))';
		T_sq_U_con = (...
			U*sqrt(diag(sigma_sq_U.^(- 1)))...
		)*(Q*pinv(Q'*Q))';
		% Contributions to Hotelling's T^2 for T: direct path
		T_sq_T_con_dir = (...
			T*sqrt(diag(sigma_sq_T.^(- 1)))...
		)*W_ast';
		% Contributions to SPE on Y
		SPE_Y_con = R_fd;
		% Contributions to SRE on X and Y
		SRE_X_con = E_fd;
		SRE_Y_con = F_fd;
	case 'absolute'
		% Upper and lower confidence limits
		lim_con = lim;
		% Contributions to Hotelling's T^2 for T and U: reconstruction path
		T_sq_T_con = (...
			T*sqrt(diag(sigma_sq_T.^(- 1)))...
		)*(P*pinv(P'*P))'.*X;
		T_sq_U_con = (...
			U*sqrt(diag(sigma_sq_U.^(- 1)))...
		)*(Q*pinv(Q'*Q))'.*Y;
		% Contributions to Hotelling's T^2 for T: direct path
		T_sq_T_con_dir = (...
			T*sqrt(diag(sigma_sq_T.^(- 1)))...
		)*W_ast'.*X;
		% Contributions to SPE on Y
		SPE_Y_con = R_fd.^2;
		% Contributions to SRE on X and Y
		SRE_X_con = E_fd.^2;
		SRE_Y_con = F_fd.^2;
end

% Loop over quality attributes
for v = 1:V_Y
	% Variable importance in projection scores
	VIP(:, v) = sqrt(...
		V_X...
		*sum(...
			max(R_sq_Y(:, v), 0)'...
			.*((W_ast/diag(sqrt(diag(W_ast'*W_ast)))).^2)...
		, 2)...
		/CR_sq_Y(end, v)...
	);
end

% Leverages of samples
obs_lev = diag((T/(T'*T))*T' + 1/N);
% Studentized residuals
MSE = sum(R.^2)./(N - dof);
R_stud = R./sqrt(MSE.*(1 - obs_lev));
% Cook's distance
Cook_D = ((R.^2)./(dof*MSE)).*(obs_lev./(1 - obs_lev));

% Confidence limit of leverage of samples
lim_obs_lev = 3*dof/N;
% Confidence limit of studentized residuals
lim_R_stud = 2;
% Confidence limit of Cook's distance
lim_Cook_D = finv(0.5, dof, N - dof);

% Assigments to the model structure
model.diagnostics.T_sq_T = T_sq_T;
model.diagnostics.T_sq_U = T_sq_U;
model.diagnostics.SPE_Y = SPE_Y;
model.diagnostics.SRE_X = SRE_X;
model.diagnostics.SRE_Y = SRE_Y;
model.diagnostics.T_sq_T_con = T_sq_T_con;
model.diagnostics.T_sq_U_con = T_sq_U_con;
model.diagnostics.T_sq_T_con_dir = T_sq_T_con_dir;
model.diagnostics.SPE_Y_con = SPE_Y_con;
model.diagnostics.SRE_X_con = SRE_X_con;
model.diagnostics.SRE_Y_con = SRE_Y_con;
model.diagnostics.VIP = VIP;
model.diagnostics.obs_lev = obs_lev;
model.diagnostics.R_stud = R_stud;
model.diagnostics.Cook_D = Cook_D;
model.diagnostics.lim_obs_lev = lim_obs_lev;
model.diagnostics.lim_R_stud = lim_R_stud;
model.diagnostics.lim_Cook_D = lim_Cook_D;

%% Estimation of confidence limits

% Choice of the method for T_sq_T and T_sq_U confidence limits
switch Tsq_lim_method
	case 'chisq'
		% chi^2 distribution approach
		lim_T_sq_T = chi_sq_limit(T_sq_T, 1 - lim);
		lim_T_sq_U = chi_sq_limit(T_sq_U, 1 - lim);
	case 'F'
		% F distribution approach
		lim_T_sq_T = F_limit(T, 1 - lim, dof);
		lim_T_sq_U = lim_T_sq_T; % same equation as the one for T
	case 'KDE'
		% Kernel density estimation approach
		lim_T_sq_T = kde_limit(T_sq_T, 1 - lim,...
			'Test_type', 'one_tail',...
			'Bandwidth', 'scott'...
		);
		lim_T_sq_U = kde_limit(T_sq_U, 1 - lim,...
			'Test_type', 'one_tail',...
			'Bandwidth', 'scott'...
		);
end

% Choice of the method for SPE_Y, SRE_X and SRE_Y confidence limits
switch SqE_lim_method
	case 'chisq'
		% chi^2 distribution approach
		lim_SPE_Y = chi_sq_limit(SPE_Y, 1 - lim);
		lim_SRE_X = chi_sq_limit(SRE_X, 1 - lim);
		lim_SRE_Y = chi_sq_limit(SRE_Y, 1 - lim);
	case 'jack_mud'
		% Jackson-Mudholkar approach
		lim_SPE_Y = jackson_mudholkar_limit(R_fd, 1 - lim);
		lim_SRE_X = jackson_mudholkar_limit(E_fd, 1 - lim);
		lim_SRE_Y = jackson_mudholkar_limit(F_fd, 1 - lim);
	case 'KDE'
		% Kernel density estimation approach
		lim_SPE_Y = kde_limit(SPE_Y, 1 - lim,...
			'Test_type', 'one_tail',...
			'Bandwidth', 'scott'...
		);
		lim_SRE_X = kde_limit(SRE_X, 1 - lim,...
			'Test_type', 'one_tail',...
			'Bandwidth', 'scott'...
		);
		lim_SRE_Y = kde_limit(SRE_Y, 1 - lim,...
			'Test_type', 'one_tail',...
			'Bandwidth', 'scott'...
		);
end

% Choice of the method for confidence limits of contributions
switch con_lim_method
	case 'norm'
		% Confindence limits for contribution to Hotelling's T^2 of T and U: reconstruction path
		lim_T_sq_T_con = norminv(lim_con, zeros(1, V_X), std(T_sq_T_con));
		lim_T_sq_U_con = norminv(lim_con, zeros(1, V_Y), std(T_sq_U_con));
		% Confindence limits for contributions to Hotelling's T^2 for T: direct path
		lim_T_sq_T_con_dir = norminv(lim_con, zeros(1, V_X), std(T_sq_T_con_dir));
		% Confindence limits for contributions to SPE on Y
		lim_SPE_Y_con = norminv(lim_con, zeros(1, V_Y), std(SPE_Y_con));
		% Confindence limits for contributions to SRE on X and Y
		lim_SRE_X_con = norminv(lim_con, zeros(1, V_X), std(SRE_X_con));
		lim_SRE_Y_con = norminv(lim_con, zeros(1, V_Y), std(SPE_Y_con));
	case 't'
		% Degrees of freedom of the distribution
		DOF = N - dof;
		% Critical t-value
		t_cl = tinv(lim_con, DOF);
		% Confindence limits for contribution to Hotelling's T^2 of T and U: reconstruction path
		lim_T_sq_T_con = sqrt(diag(T_sq_T_con'*T_sq_T_con)/DOF)'*t_cl;
		lim_T_sq_U_con = sqrt(diag(T_sq_U_con'*T_sq_U_con)/DOF)'*t_cl;
		% Confindence limits for contributions to Hotelling's T^2 for T: direct path
		lim_T_sq_T_con_dir = sqrt(diag(T_sq_T_con_dir'*T_sq_T_con_dir)/DOF)'*t_cl;
		% Confindence limits for contributions to SPE on Y
		lim_SPE_Y_con = sqrt(diag(SPE_Y_con'*SPE_Y_con)/DOF)'*t_cl;
		% Confindence limits for contributions to SRE on X and Y
		lim_SRE_X_con = sqrt(diag(SRE_X_con'*SRE_X_con)/DOF)'*t_cl;
		lim_SRE_Y_con = sqrt(diag(SRE_Y_con'*SRE_Y_con)/DOF)'*t_cl;
% 	case 'quantile'
% 		% Confindence limits for contribution to Hotelling's T^2 of T and U: reconstruction path
% 		lim_T_sq_T_con = quantile(T_sq_T_con, lim_con);
% 		lim_T_sq_U_con = quantile(T_sq_U_con, lim_con);
% 		% Confindence limits for contributions to Hotelling's T^2 for T: direct path
% 		lim_T_sq_T_con_dir = quantile(T_sq_T_con_dir, lim_con);
% 		% Confindence limits for contributions to SPE on Y
% 		lim_SPE_Y_con = quantile(SPE_Y_con, lim_con);
% 		% Confindence limits for contributions to SRE on X and Y
% 		lim_SRE_X_con = quantile(SRE_X_con, lim_con);
% 		lim_SRE_Y_con = quantile(SRE_Y_con, lim_con);
end

% Choice of the semiaxes for confidence ellipse
switch l_kind
	case 'two'
		% Issue a warning is chi^2 method is used for confidence limits
		if ~strcmp(Tsq_lim_method, 'F')
			warning(['Confidence ellipse requested as based on two latent '...
				'variables only, possible only using F distribution-based '...
				'confidence limits for T_sq_T and T_sq_U: limits will be '...
				'inconsistent'])
		end
		lim_T_sq_T_for_plots = F_limit(T, 1 - lim, 2);
		lim_T_sq_U_for_plots = lim_T_sq_T_for_plots;
	case 'full'
		lim_T_sq_T_for_plots = lim_T_sq_T;
		lim_T_sq_U_for_plots = lim_T_sq_U;
end
% Confidence ellipsoid semiaxes of T and U
l_T = sqrt(sigma_sq_T*lim_T_sq_T_for_plots);
l_U = sqrt(sigma_sq_U*lim_T_sq_U_for_plots);

% Assigments to the model structure
model.estimates.lim = lim;
model.estimates.dof = dof;
model.estimates.lim_T_sq_T = lim_T_sq_T;
model.estimates.lim_T_sq_U = lim_T_sq_U;
model.estimates.lim_SPE_Y = lim_SPE_Y;
model.estimates.lim_SRE_X = lim_SRE_X;
model.estimates.lim_SRE_Y = lim_SRE_Y;
model.estimates.lim_T_sq_T_con = lim_T_sq_T_con;
model.estimates.lim_T_sq_U_con = lim_T_sq_U_con;
model.estimates.lim_T_sq_T_con_dir = lim_T_sq_T_con_dir;
model.estimates.lim_SPE_Y_con = lim_SPE_Y_con;
model.estimates.lim_SRE_X_con = lim_SRE_X_con;
model.estimates.lim_SRE_Y_con = lim_SRE_Y_con;
model.estimates.l_T = l_T;
model.estimates.l_U = l_U;

%% Output assignments

model_out = model;

end