function CV_out = cross_validate_pls (X_in, Y_in, A_in, method_in, varargin)
% CROSS_VALIDATE_PLS Cross validation of PLS model with population analysis
%
% Syntax:
%	CV = cross_validate_pls (X, Y, A, 'leave_one_out')
%	CV = cross_validate_pls (X, Y, A, 'continuous_blocks', G)
%	CV = cross_validate_pls (X, Y, A, 'venetian_blind', G)
%	CV = cross_validate_pls (X, Y, A, 'venetian_blind', G, band_thick)
%	CV = cross_validate_pls (X, Y, A, 'random_subsets', G, N_rep)
%	CV = cross_validate_pls (_, 'Key', value)
%
% Inputs:
%	X:			Structure of arrays of predictor variables
%	Y:			Structure of arrays of response variables
%	A:			Number of latent variables to assess in cross-validation
%	method:		Method to be used for cross-validation (grouping of observations)
%	G:			Number of groups to be generated (only for continuous blocks,
%				venetian blind and random subsets methods)
%	band_thick:	Thickenss of a single band of the blind (only for venetian blind
%				method, optional, default = 1)
%	N_rep:		Number of iterations (only for random subsets method)
%
% Outputs:
%	CV:	Structure containing results of cross-validation
%
% Keys and Values:
%	'PrintInfo': [true | false]
%		Print progress of cross-validatioon to the command window
%		(default = false)
%	'Aggregate': [true | false]
%		Whether to aggregate the MSE used for computing of the suggested number
%		of latent variables for all Y variables; by default, a number of latent
%		variables is suggestd for each Y vairbale separately (false), while
%		aggregation finds a comprimise suggested number of latent variables.
%	'ParFac': [pfac]
%		Parsimony factor for standard error rule used to compute the suggested
%		number of latent variables: setting pfac = 0 suggests the number of
%		latent variables which minimise the MSE, while setting pfac = 1 (default)
%		uses the 1-standard error rule; any other value is in principle
%		admissible.
%	'Preprocessing': ['none' | 'mean_centre' | 'autoscale']
%		Preprocessing method to be applied to the input data structures, can be:
%		no scaling at all ('none'); mean-centring only ('mean_centre');
%		mean-centring and scaling to unit variance ('autoscale', default); note
%		that the algorithm implemented by this function requires mean_centred
%		data matrices, therefore if the option 'none' is used the user is in
%		charge of providing mean_centred matrices, otherwise an error is issued
%		(default = 'autoscale')
%	'Algorithm': ['simpls' | 'nipals']
%		Algorithm to be used for PLS calibration, can be either
%		statistically-inspired modification to partial least squares (SIMPLS) or
%		non-linear iteratrive partial least squared (NIPaLS) (default = 'nipals')
%	'Tol': [tol]
%		Tolerance for convergence as 2-norm on relative scores variation
%		(NIPaLS algorithm only, default = 1e-11)
%	'MaxIter': [max_iter]
%		Maximum number of iterations allowed (NIPaLS algorithm only,
%		default = 150)
%	'Normalise': ['standard' | 'loadings']
%		Normalisation scheme to be adopted, either standard normalisation
%		according to the algorithm of choice or PCA-coherent normalisation,
%		meaning both the X loadings and Y loadings are normalised (default =
%		'standard'); see functions pls_calibration_simpls and
%		pls_calibration_nipals for more information on normalisation schemes
%	'ConLimCV': [lim_CV]
%		Confidence limit for statistics on model parameters for estimating
%		confidence limits on model parameters in cross-validation (default = 0.95)
%	'DOFMethodCV: ['naive' | 'cv_based' | 'krylov']
%		Method for computing degrees of freedom of the model for estimation of
%		confidence limits on model parameters and on estimates of control limits,
%		can be either based on the number of latent variables ('naive',
%		dof = A, default), based on leave-one-out cross-validation (cv_based,
%		dof estimated according to Van Der Voet, 1999) or based on Krylov
%		sequence ('krylov', dof estimated according to Kramer, 2012)
%	'LimMethodCV': ['norm' | 't']
%		Method for computing the confidence limits of statistics on model
%		parameters in cross-validation, can be based on a normal distribution or
%		on a t distribution (default = 'norm')
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
%	The 'ErrBasedOn' key is not provided in this function as PRESS, RMSECV,
%	CV_bias and SECV will be reported as scaled even in the unscaled errors are
%	requested as they are meant to assess the performance of the model and the
%	comparison of errors on different variables is possible only if they are all
%	scaled

%{
cross_validate_pls.m
Version: 1.0.1
Date: 2023-12-20
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

cross_validate_pls.m
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

X = X_in;
Y = Y_in;
A = A_in;
method = method_in;

%% Initial checks

% Check if the data array is unfolded
if size(X, 3) ~= 1 || size(Y, 3) ~= 1
	error('The data arrays must be unfolded for PLS model calibration')
end
% Check if the requested number of LVs is feasible
if A > min(size(X))
	error(['The number of latent variables cannot exceed min(size(X)) = ...'
		num2str(min(size(X)))])
end
% Check if there are missing values
if sum(ismissing(X), 'all') ~= 0 || sum(ismissing(Y), 'all') ~= 0
	error('Missing values found: PLS cannot be calibrated')
end

% Number of observations and number of variables
[N, V_X] = size(X);
V_Y = size(Y, 2);

%% Optional arguments development

% Optionals initialization
G_obs = [];
band_thick_obs = [];
N_rep = 1;
printinfo = false;
aggr = false;
pfac = 1;
preprocess = 'autoscale';
alg = 'nipals';
tol = 1e-11;
max_iter = 150;
norm_sch = 'standard';
err_on = 'scaled';
lim_CV = 0.95;
dof_method_CV = 'naive';
CV_lim_method = 'norm';

obs_names = cellstr([repmat('O', N, 1) char(pad(replace(string(num2str((1:N)')),  ' ', ''), length(num2str(N)), 'left', '0'))]);
X_var_names = cellstr([repmat('X', V_X, 1) char(pad(replace(string(num2str((1:V_X)')),  ' ', ''), length(num2str(V_X)), 'left', '0'))]);
Y_var_names = cellstr([repmat('Y', V_Y, 1) char(pad(replace(string(num2str((1:V_Y)')),  ' ', ''), length(num2str(V_Y)), 'left', '0'))]);

% Development cycle
if ~isempty(varargin)
	if isnumeric(varargin{1})
		G_obs = varargin{1};
		if length(varargin) > 1 && isnumeric(varargin{2})
			switch method
				case 'venetian_blind'
					band_thick_obs = varargin{2};
				case 'random_subsets'
					N_rep = varargin{2};
			end
			varargin = varargin(3:end);
		else
			varargin = varargin(2:end);
		end
	end
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
			case 'PrintInfo'
				printinfo = varargin{i + 1};
			case 'Aggregate'
				aggr = varargin{i + 1};
			case 'ParFac'
				pfac = varargin{i + 1};
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
			case 'ConLimCV'
				lim_CV = varargin{i + 1};
			case 'DOFMethodCV'
				dof_method_CV = varargin{i + 1};
				if strcmp(dof_method_CV, 'naive') || strcmp(dof_method_CV, 'cv_based') || strcmp(dof_method_CV, 'krylov')
				else
					error(['Supported methods for degrees of freedom '...
						'estimation are naive, cv_based and krylov'])
				end
			case 'LimMethodCV'
				CV_lim_method = varargin{i + 1};
				if strcmp(CV_lim_method, 'norm') || strcmp(CV_lim_method, 't')
				else
					error(['Supported methods for limits are normal '...
						'distribution method or t distribution method'])
				end
			case 'DiagBasedOn'
				diag_on = varargin{i + 1};
				if strcmp(diag_on, 'scaled') || strcmp(diag_on, 'unscaled')
				else
					error('Undefined key for DiagBasedOn')
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
if strcmp(dof_method_CV, 'cv_based') && V_Y > 1
	warning(['cv_based DOF estimation is meant for univariate Y: a variation '...
		'for multivariate Y will be used. However, for multivariate Y it is '...
		'suggested to use naive DOF estimation'])
end

%% Initialise results structure

CV = struct;
	% Results
	CV.results = struct;
		CV.results.A_sug = [];
		CV.results.R_sq_Y = [];
		CV.results.PRESS = [];
		CV.results.RMSECV = [];
		CV.results.CV_bias = [];
		CV.results.SECV = [];
	% Statistics on model parameters
	CV.statistics = struct;
		CV.statistics.lim = [];
		CV.statistics.dof = [];
		CV.statistics.R_sq_Y = [];
		CV.statistics.PRESS = [];
		CV.statistics.RMSECV = [];
		CV.statistics.CV_bias = [];
		CV.statistics.SECV = [];
		CV.statistics.R = [];
		CV.statistics.W_ast = [];
		CV.statistics.P = [];
		CV.statistics.Q = [];
		CV.statistics.B = [];
		CV.statistics.VIP = [];
		CV.statistics.MSE_split = [];
		% add fields for diagnostics
	% Dimensions of the model entities
	CV.dimensions = struct;
		CV.dimensions.N = [];
		CV.dimensions.V_X = [];
		CV.dimensions.V_Y = [];
		CV.dimensions.A = [];
		CV.dimensions.N_rep = [];
	% Infos on cross-validation
	CV.info = struct;
		CV.info.CV_method = '';
		CV.info.grouping = [];
		CV.info.N_rep = [];
		CV.info.pfac = [];
		CV.info.preprocessing = '';
		CV.info.normalisation = '';
		CV.info.error_based_on = '';
		CV.info.dof_method_CV = '';
		CV.info.obs_names = {};
		CV.info.X_var_names = {};
		CV.info.Y_var_names = {};

% Dimesions and labels assignments
CV.dimensions.N = N;
CV.dimensions.V_X = V_X;
CV.dimensions.V_Y = V_Y;
CV.dimensions.A = A;
CV.dimensions.N_rep = N_rep;
CV.info.pfac = pfac;
CV.info.preprocessing = preprocess;
CV.info.normalisation = norm_sch;
CV.info.error_based_on = err_on;
CV.info.dof_method_CV = dof_method_CV;
CV.info.obs_names = obs_names;
CV.info.X_var_names = X_var_names;
CV.info.Y_var_names = Y_var_names;

%% DOF of the model

if printinfo
	fprintf('Estimating DOF of the model\n')
	tic
end

% Degrees of freedom of the model
[dof_CV, ~] = pls_dof(X, Y, A,...
	'Method', dof_method_CV,...
	'Preprocessing', preprocess,...
	'Algorithm', alg,...
	'Tol', tol,...
	'MaxIter', max_iter,...
	'Normalise', norm_sch...
);

% Compress DOF for multivariate Y
dof_CV = sum(dof_CV, 2)/V_Y;

%% Iteration loop

% Pre-allocations
R_sq_CV_reps = zeros(A, V_Y, N_rep);
PRESS_reps = zeros(A, V_Y, N_rep);
CV_bias_reps = zeros(A, V_Y, N_rep);
SECV_reps = zeros(A, V_Y, N_rep);

R_CV_pop_reps = zeros(N, V_Y, A, N_rep);

if strcmp(method, 'leave_one_out')
	W_ast_pop_reps = zeros(V_X, A, N, N_rep);
	P_pop_reps = zeros(V_X, A, N, N_rep);
	Q_pop_reps = zeros(V_Y, A, N, N_rep);
	B_pop_reps = zeros(V_X, V_Y, A, N, N_rep);
	VIP_pop_reps = zeros(V_X, V_Y, A, N, N_rep);

	MSE_split_reps = zeros(A, V_Y, N, N_rep);
else
	W_ast_pop_reps = zeros(V_X, A, G_obs, N_rep);
	P_pop_reps = zeros(V_X, A, G_obs, N_rep);
	Q_pop_reps = zeros(V_Y, A, G_obs, N_rep);
	B_pop_reps = zeros(V_X, V_Y, A, G_obs, N_rep);
	VIP_pop_reps = zeros(V_X, V_Y, A, G_obs, N_rep);

	MSE_split_reps = zeros(A, V_Y, G_obs, N_rep);
end

grouping_reps = zeros(N, N_rep);

% Print info
if printinfo && strcmp(method, 'random_subsets')
	fprintf('Pre-allocations done, starting cross-validation\n')
	tic
end

% iteration loop
for i = 1:N_rep
	
	% Print info
	if printinfo && strcmp(method, 'random_subsets')
		fprintf('CV iteration: %.3g/%.3g\n', i, N_rep)
	end

	%% Grouping of samples

	% Grouping of observations
	grouping_obs = cross_validation_grouping(method, N, G_obs, band_thick_obs);
	% Ensure that G_obs is not empty
	G_obs = max(grouping_obs);
	
	% Save grouping outer variable
	grouping_reps(:, i) = grouping_obs;

	%% Utilities and preallocations

	% Pre-allocations
	R_CV = zeros(N, V_Y, A);
	E_CV = zeros(N, V_X, A);
	F_CV = zeros(N, V_Y, A);
	
	R_CV_uns = zeros(N, V_Y, A);
	
	W_ast_pop = zeros(V_X, A, G_obs);
	P_pop = zeros(V_X, A, G_obs);
	Q_pop = zeros(V_Y, A, G_obs);
	B_pop = zeros(V_X, V_Y, A, G_obs);
	VIP_pop = zeros(V_X, V_Y, A, G_obs);

	R_sq_aux = zeros(A, V_Y);
	
	% Print info
	if printinfo && ~strcmp(method, 'random_subsets')
		fprintf('Pre-allocations done, starting cross-validation\n')
		tic
	end

	%% Cross-validation loop

	for g = 1:G_obs
		
		% Print progress of CV
		if printinfo && ~strcmp(method, 'random_subsets')
			fprintf('CV group: %.3g/%.3g\n', g, G_obs)
		end
		
		% Validation group index
		gvi = g == grouping_obs;
		% Calibration dataset
		X_cal = X(~gvi, :);
		Y_cal = Y(~gvi, :);
		% Validation sample
		X_val = X(gvi, :);
		Y_val = Y(gvi, :);
		
		% Scale calibration data
		switch preprocess
			case 'none'
				mu_X = zeros(1, V_X);
				sigma_X = ones(1, V_X);
				X_cal_s = X_cal;
				mu_Y = zeros(1, V_Y);
				sigma_Y = ones(1, V_Y);
				Y_cal_s = Y_cal;
			case 'mean_centre'
				mu_X = mean(X_cal);
				sigma_X = ones(1, V_X);
				X_cal_s = X_cal - mu_X;
				mu_Y = mean(Y_cal);
				sigma_Y = ones(1, V_Y);
				Y_cal_s = Y_cal - mu_Y;
			case 'autoscale'
				[X_cal_s, mu_X, sigma_X] = autoscale(X_cal);
				[Y_cal_s, mu_Y, sigma_Y] = autoscale(Y_cal);
		end
		% Scale validation data
		X_val_s = scale_by(X_val, mu_X, sigma_X);
		Y_val_s = scale_by(Y_val, mu_Y, sigma_Y);
		
		% Build the model
		switch alg
			case 'simpls'
				[~, W_ast, P, T, Q, ~, b] = pls_calibration_simpls (X_cal_s, Y_cal_s, A, 'Normalise', norm_sch);
			case 'nipals'
				[~, W_ast, P, T, Q, ~, b] = pls_calibration_nipals(X_cal_s, Y_cal_s, A, 'Tol', tol, 'MaxIter', max_iter, 'Normalise', norm_sch);
		end
		
		% Total variance for each Y for R_sq
		ssq_Y_pv = diag(Y_cal_s'*Y_cal_s)';
		% Compute scores applying the model
		T_val = X_val_s*W_ast;
		U_val = Y_val_s*Q*pinv(Q'*Q);
		% Loop over LVs
		for a = 1:A
			% Residuals in cross-validation
			R_CV(gvi, :, a) = Y_val_s - T_val(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
			% Reconstruction errors in cross-validation
			E_CV(gvi, :, a) = X_val_s - T_val(:, 1:a)*P(:, 1:a)';
			F_CV(gvi, :, a) = Y_val_s - U_val(:, 1:a)*Q(:, 1:a)';
			% Determination coefficient on each response variable in calibration
			R_sq_aux(a, :) = 1 - sum((Y_cal_s - T(:, a)*b(a)*Q(:, a)').^2, 1)./ssq_Y_pv;
			% Unscaled cross-validation residuals for determination coefficients
			R_CV_uns(gvi, :, a) = R_CV(gvi, :, a)*diag(sigma_Y);
		end
		% Cumulative determination coefficient on each response variable in calibration
		CR_sq_aux = cumsum(R_sq_aux);

		% Save samples of populations of parameters
		W_ast_pop(:, :, g) = W_ast;
		P_pop(:, :, g) = P;
		Q_pop(:, :, g) = Q;
		
		% Loop over LVs
		for a = 1:A
			% PLS regressione coefficients
			B_pop(:, :, a, g) = W_ast(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
			% Loop over quality attributes
			for v = 1:V_Y
				% Variable importance in projection scores
				VIP_pop(:, v, a, g) = sqrt(...
					V_X...
					*sum(...
						max(R_sq_aux(1:a, v), 0)'...
						.*((W_ast(:, 1:a)/diag(sqrt(diag(W_ast(:, 1:a)'*W_ast(:, 1:a))))).^2)...
					, 2)...
					/CR_sq_aux(a, v)...
				);
			end
		end
		% Assignment to MSE population
		MSE_split_reps(:, :, g, i) = permute(sum((R_CV(gvi, :, :)).^2, 1), [3, 2, 1])/sum(gvi);
	end
	
	% Determination coefficient in cross-validation
	R_sq_CV = 1 - permute(sum((R_CV_uns).^2, 1), [3, 2, 1])./sum((Y - mean(Y)).^2);

	% Prediction error sum of squares update
	PRESS = permute(sum(R_CV.^2, 1), [3, 2, 1]);
	% Cross-validation error bias
	CV_bias = permute(sum(R_CV, 1)/N, [3, 2, 1]);
	% Stadard error in cross-validation
	SECV = permute(sqrt(sum((R_CV - sum(R_CV, 1)/N).^2, 1)/(N - 1)), [3, 2, 1]);

	% Save samples of populations of parameters in outer variables
	R_sq_CV_reps(:, :, i) = R_sq_CV;
	PRESS_reps(:, :, i) = PRESS;
	CV_bias_reps(:, :, i) = CV_bias;
	SECV_reps(:, :, i) = SECV;
	R_CV_pop_reps(:, :, :, i) = R_CV;
	W_ast_pop_reps(:, :, :, i) = W_ast_pop;
	P_pop_reps(:, :, :, i) = P_pop;
	Q_pop_reps(:, :, :, i) = Q_pop;
	B_pop_reps(:, :, :, :, i) = B_pop;
	VIP_pop_reps(:, :, :, :, i) = VIP_pop;
end

% Root mean squared error in cross-validation
RMSECV_reps = sqrt(PRESS_reps/N);

% Assigments to the model structure
CV.info.CV_method = method;
CV.info.grouping = grouping_reps;
CV.info.N_rep = N_rep;

%% Aggregate over repetitions

% Determination coefficient in cross-validation
R_sq_CV = sum(R_sq_CV_reps, 3)/N_rep;
% Root mean squared error in cross-validation
PRESS = sum(PRESS_reps, 3)/N_rep;
% Root mean squared error in cross-validation
RMSECV = sum(RMSECV_reps, 3)/N_rep;
% Cross-validation error bias
CV_bias =  sum(CV_bias_reps, 3)/N_rep;
% Stadard error of prediction in cross-validation
SECV =  sum(SECV_reps, 3)/N_rep;

% Assigments to the model structure
CV.results.R_sq_Y = R_sq_CV;
CV.results.PRESS = PRESS;
CV.results.RMSECV = RMSECV;
CV.results.CV_bias = CV_bias;
CV.results.SECV = SECV;

%% Analysis of populations

% Print info
if printinfo
	fprintf('\nCross-validation done, estimating variability of parameters\n\n')
end

% Two-tail confidence limits
u_lim_CV = lim_CV + (1 - lim_CV)/2;

% Assigments to the model structure: poulations
CV.statistics.R_sq_Y.pop = R_sq_CV_reps;
CV.statistics.PRESS.pop = PRESS_reps;
CV.statistics.RMSECV.pop = RMSECV_reps;
CV.statistics.CV_bias.pop = CV_bias_reps;
CV.statistics.SECV.pop = SECV_reps;
CV.statistics.R.pop = R_CV_pop_reps;
CV.statistics.W_ast.pop = W_ast_pop_reps;
CV.statistics.P.pop = P_pop_reps;
CV.statistics.Q.pop = Q_pop_reps;
CV.statistics.B.pop = B_pop_reps;
CV.statistics.VIP.pop = VIP_pop_reps;
% Assigments to the model structure: means
CV.statistics.R_sq_Y.mean = mean(R_sq_CV_reps, 3);
CV.statistics.PRESS.mean = mean(PRESS_reps, 3);
CV.statistics.RMSECV.mean = mean(RMSECV_reps, 3);
CV.statistics.CV_bias.mean = mean(CV_bias_reps, 3);
CV.statistics.SECV.mean = mean(SECV_reps, 3);
CV.statistics.R.mean = mean(R_CV_pop_reps, 4);
CV.statistics.W_ast.mean = mean(W_ast_pop_reps, [3, 4]);
CV.statistics.P.mean = mean(P_pop_reps, [3, 4]);
CV.statistics.Q.mean = mean(Q_pop_reps, [3, 4]);
CV.statistics.B.mean = mean(B_pop_reps, [4, 5]);
CV.statistics.VIP.mean = mean(VIP_pop_reps, [4, 5]);
% Assigments to the model structure: standard deviations
CV.statistics.R_sq_Y.std = std(R_sq_CV_reps, 0, 3);
CV.statistics.PRESS.std = std(PRESS_reps, 0, 3);
CV.statistics.RMSECV.std = std(RMSECV_reps, 0, 3);
CV.statistics.CV_bias.std = std(CV_bias_reps, 0, 3);
CV.statistics.SECV.std = std(SECV_reps, 0, 3);
CV.statistics.R.std = std(R_CV_pop_reps, 0, 4);
CV.statistics.W_ast.std = std(W_ast_pop_reps, 0, [3, 4]);
CV.statistics.P.std = std(P_pop_reps, 0, [3, 4]);
CV.statistics.Q.std = std(Q_pop_reps, 0, [3, 4]);
CV.statistics.B.std = std(B_pop_reps, 0, [4, 5]);
CV.statistics.VIP.std = std(VIP_pop_reps, 0, [4, 5]);
% Assigments to the model structure: limits around a null value
switch CV_lim_method
	case 'norm'
		CV.statistics.R_sq_Y.cl = norminv(u_lim_CV, zeros(A, V_Y), CV.statistics.R_sq_Y.std);
		CV.statistics.PRESS.cl = norminv(u_lim_CV, zeros(A, V_Y), CV.statistics.PRESS.std);
		CV.statistics.RMSECV.cl = norminv(u_lim_CV, zeros(A, V_Y), CV.statistics.RMSECV.std);
		CV.statistics.CV_bias.cl = norminv(u_lim_CV, zeros(A, V_Y), CV.statistics.CV_bias.std);
		CV.statistics.SECV.cl = norminv(u_lim_CV, zeros(A, V_Y), CV.statistics.SECV.std);
		CV.statistics.R.cl = norminv(u_lim_CV, zeros(N, V_Y, A), CV.statistics.R.std);
		CV.statistics.W_ast.cl = norminv(u_lim_CV, zeros(V_X, A), CV.statistics.W_ast.std);
		CV.statistics.P.cl = norminv(u_lim_CV, zeros(V_X, A), CV.statistics.P.std);
		CV.statistics.Q.cl = norminv(u_lim_CV, zeros(V_Y, A), CV.statistics.Q.std);
		CV.statistics.B.cl = norminv(u_lim_CV, zeros(V_X, V_Y, A), CV.statistics.B.std);
		CV.statistics.VIP.cl = norminv(u_lim_CV, zeros(V_X, V_Y, A), CV.statistics.VIP.std);
	case 't'
		% Pre-allocations
		CV.statistics.R_sq_Y.cl = zeros(A, V_Y);
		CV.statistics.PRESS.cl = zeros(A, V_Y);
		CV.statistics.RMSECV.cl = zeros(A, V_Y);
		CV.statistics.CV_bias.cl = zeros(A, V_Y);
		CV.statistics.SECV.cl = zeros(A, V_Y);
		CV.statistics.R.cl = zeros(N, V_Y, A);
		CV.statistics.W_ast.cl = zeros(V_X, A);
		CV.statistics.P.cl = zeros(V_X, A);
		CV.statistics.Q.cl = zeros(V_Y, A);
		CV.statistics.B.cl = zeros(V_X, V_Y, A);
		CV.statistics.VIP.cl = zeros(V_X, V_Y, A);
		% Loop over LVs
		for a = 1:A
			DOF = max(N_rep - dof_CV(a), 0);
			t_cl = tinv(u_lim_CV, DOF);
			CV.statistics.R_sq_Y.cl(a, :) = sqrt(sum((CV.statistics.R_sq_Y.pop(a, :, :) - CV.statistics.R_sq_Y.mean(a, :)).^2, 3)/DOF)*t_cl;
			CV.statistics.PRESS.cl(a, :) = sqrt(sum((CV.statistics.PRESS.pop(a, :, :) - CV.statistics.PRESS.mean(a, :)).^2, 3)/DOF)*t_cl;
			CV.statistics.RMSECV.cl(a, :) = sqrt(sum((CV.statistics.RMSECV.pop(a, :, :) - CV.statistics.RMSECV.mean(a, :)).^2, 3)/DOF)*t_cl;
			CV.statistics.CV_bias.cl(a, :) = sqrt(sum((CV.statistics.CV_bias.pop(a, :, :) - CV.statistics.CV_bias.mean(a, :)).^2, 3)/DOF)*t_cl;
			CV.statistics.SECV.cl(a, :) = sqrt(sum((CV.statistics.SECV.pop(a, :, :) - CV.statistics.SECV.mean(a, :)).^2, 3)/DOF)*t_cl;
			CV.statistics.R.cl(:, :, a) = sqrt(sum((CV.statistics.R.pop(:, :, a, :) - CV.statistics.R.mean(:, :, a)).^2, 4)/DOF)*t_cl;
			DOF = max(G_obs*N_rep - dof_CV(a), 0);
			t_cl = tinv(u_lim_CV, DOF);
			CV.statistics.W_ast.cl(:, a) = sqrt(sum((CV.statistics.W_ast.pop(:, a, :, :) - CV.statistics.W_ast.mean(:, a)).^2, [3, 4])/DOF)*t_cl;
			CV.statistics.P.cl(:, a) = sqrt(sum((CV.statistics.P.pop(:, a, :, :) - CV.statistics.P.mean(:, a)).^2, [3, 4])/DOF)*t_cl;
			CV.statistics.Q.cl(:, a) = sqrt(sum((CV.statistics.Q.pop(:, a, :, :) - CV.statistics.Q.mean(:, a)).^2, [3, 4])/DOF)*t_cl;
			CV.statistics.B.cl(:, :, a) = sqrt(sum((CV.statistics.B.pop(:, :, a, :, :) - CV.statistics.B.mean(:, :, a)).^2, [4, 5])/DOF)*t_cl;
			CV.statistics.VIP.cl(:, :, a) = sqrt(sum((CV.statistics.VIP.pop(:, :, a, :, :) - CV.statistics.VIP.mean(:, :, a)).^2, [4, 5])/DOF)*t_cl;
		end
% 	case 'quantile'
% 		CV.statistics.R_sq_Y.lcl = quantile(CV.statistics.R_sq_Y.pop, l_lim_CV, 3);
% 		CV.statistics.PRESS.lcl = quantile(CV.statistics.PRESS.pop, l_lim_CV, 3);
% 		CV.statistics.RMSECV.lcl = quantile(CV.statistics.RMSECV.pop, l_lim_CV, 3);
% 		CV.statistics.CV_bias.lcl = quantile(CV.statistics.CV_bias.pop, l_lim_CV, 3);
% 		CV.statistics.SECV.lcl = quantile(CV.statistics.SECV.pop, l_lim_CV, 3);
% 		CV.statistics.R.lcl = quantile(CV.statistics.R.pop, l_lim_CV, 4);
% 		
% 		CV.statistics.W_ast.lcl = quantile(CV.statistics.W_ast.pop, l_lim_CV, [3, 4]);
% 		CV.statistics.P.lcl = quantile(CV.statistics.P.pop, l_lim_CV, [3, 4]);
% 		CV.statistics.Q.lcl = quantile(CV.statistics.Q.pop, l_lim_CV, [3, 4]);
% 		CV.statistics.B.lcl = quantile(CV.statistics.B.pop, l_lim_CV, [4, 5]);
% 		CV.statistics.VIP.lcl = quantile(CV.statistics.VIP.pop, l_lim_CV, [4, 5]);
% 		CV.statistics.R_sq_Y.ucl = quantile(CV.statistics.R_sq_Y.pop, u_lim_CV, 3);
% 		CV.statistics.PRESS.ucl = quantile(CV.statistics.PRESS.pop, u_lim_CV, 3);
% 		CV.statistics.RMSECV.ucl = quantile(CV.statistics.RMSECV.pop, u_lim_CV, 3);
% 		CV.statistics.CV_bias.ucl = quantile(CV.statistics.CV_bias.pop, u_lim_CV, 3);
% 		CV.statistics.SECV.ucl = quantile(CV.statistics.SECV.pop, u_lim_CV, 3);
% 		CV.statistics.R.ucl = quantile(CV.statistics.R.pop, u_lim_CV, 4);
% 		
% 		CV.statistics.W_ast.ucl = quantile(CV.statistics.W_ast.pop, u_lim_CV, [3, 4]);
% 		CV.statistics.P.ucl = quantile(CV.statistics.P.pop, u_lim_CV, [3, 4]);
% 		CV.statistics.Q.ucl = quantile(CV.statistics.Q.pop, u_lim_CV, [3, 4]);
% 		CV.statistics.B.ucl = quantile(CV.statistics.B.pop, u_lim_CV, [4, 5]);
% 		CV.statistics.VIP.ucl = quantile(CV.statistics.VIP.pop, u_lim_CV, [4, 5]);
end

% Assigments to the model structure
CV.statistics.lim = lim_CV;
CV.statistics.dof = dof_CV;

%% Suggested number of LVs

% Selection of aggregation scheme
if aggr
	% Mean of MSE population
	mean_MSE = mean(MSE_split_reps, [2, 3, 4]);
	% Standard deviation of MSE population
	std_MSE = std(MSE_split_reps, 0, [2, 3, 4]);
	% Minimum of MSE
	min_MSE = min(mean_MSE, [], 1);
	% Pointers for the optimal values
	idx = mean_MSE <= min_MSE + (pfac*std_MSE)/sqrt(G_obs*N_rep);
	% Suggested number of LVs
	A_sug = find(idx, 1);
else
	% Mean of MSE population
	mean_MSE = mean(MSE_split_reps, [3, 4]);
	% Standard deviation of MSE population
	std_MSE = std(MSE_split_reps, 0, [3, 4]);
	% Minimum of MSE
	min_MSE = min(mean_MSE, [], 1);
	% Pointers for the optimal values
	idx = mean_MSE <= min_MSE + (pfac*std_MSE)/sqrt(G_obs*N_rep);
	% Suggested number of LVs
	A_sug = zeros(1, V_Y);
	for w = 1:V_Y
		A_sug(w) = find(idx(:, w), 1);
	end
end

% Assigments to the model structure
CV.results.A_sug = A_sug;
CV.statistics.MSE_split.pop = MSE_split_reps;
CV.statistics.MSE_split.mean = mean_MSE;
CV.statistics.MSE_split.std = std_MSE;

% Print info
if printinfo
	t_ela = toc;
	fprintf('Cross-validation completed in %.4g second\n\n', t_ela)
end

%% Output assignments

CV_out = CV;

end