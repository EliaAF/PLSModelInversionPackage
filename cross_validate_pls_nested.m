function NCV_out = cross_validate_pls_nested (X_in, Y_in, A_in, method_in, varargin)
% CROSS_VALIDATE_PLS_NESTED PLS model nested cross-validation with population
%	analysis
%
% Syntax:
%	NCV = cross_validate_pls_nested (X, Y, A, 'leave_one_out')
%	NCV = cross_validate_pls_nested (X, Y, A, 'continuous_blocks', G)
%	NCV = cross_validate_pls_nested (X, Y, A, 'venetian_blind', G)
%	NCV = cross_validate_pls_nested (X, Y, A, 'venetian_blind', G, band_thick)
%	NCV = cross_validate_pls_nested (X, Y, A, 'random_subsets', G, N_rep)
%	NCV = cross_validate_pls_nested (_, 'Key', value)
%
% Inputs:
%	X:			Structure of arrays of predictor variables
%	Y:			Structure of arrays of response variables
%	A:			Number of latent variables to assess in cross-validation outer
%				loop
%	method:		Method to be used for cross-validation outer loop (grouping of
%				observations)
%	G:			Number of groups to be generated (only for continuous blocks,
%				venetian blind and random subsets methods)
%	band_thick:	Thickenss of a single band of the blind (only for venetian blind
%				method, optional, default = 1)
%	N_rep:		Number of iterations of the outer loop (only for random subsets
%				method)
%
% Outputs:
%	NCV:	Structure containing results of nested cross-validation
%
% Keys and Values:
%	'PrintInfo': [true | false]
%		Print progress of cross-validation to the command window
%		(default = false)
%	'Aggregate': [true | false]
%		Whether to aggregate final number of latent variables for all Y
%		variables; by default, a final number of latent variables is reported for
%		each Y vairbale separately (false), while aggregation reports a
%		comprimise suggested number of latent variables. This option controls
%		also the aggregation scheme used in the cross-validation inner loop. Note
%		that this option also influences the results of the prediction error
%		estimation: without aggregation errors are assumed to be computed using
%		different models (otpimised for each one of the Y variables one at a
%		time).
%	'ParFac': [pfac]
%		Parsimony factor for standard error rule used to compute the suggested
%		number of latent variables: setting pfac = 0 suggests the number of
%		latent variables which minimise the MSE, while setting pfac = 1 (default)
%		uses the 1-standard error rule; any other value is in principle
%		admissible.
%	'InnerMethod': ['leave_one_out' | 'continuous_blocks' | 'venetian_blind' | 'random_subsets']
%		Method to be used for grouping of observations in cross-validation inner
%		loop (default = 'continuous_blocks')
%	'InnerGroups': [G_obs_inner]
%		Number of groups to be generated for grouping of observations in
%		cross-validation inner loop (defualt = G_obs - 1); this key has effect
%		only if InnerMethod is set to 'continuous_blocks', 'venetian_blind' or
%		'random_subsets'
%	'InnerBandThickness': [band_thick_inner]
%		Thickenss of a single band of the blind for grouping of observations in
%		cross-validation inner loop (defualt = 1); this key has effect only if
%		InnerMethod is set to 'venetian_blind'
%	'InnerRepetitions': [N_rep_inner]
%		Number of iterations in cross-validation inner loop; this key has effect
%		only if InnerMethod is set to 'random_subsets'
%	'InnerA': [A_inner]
%		Number of latent variables to be assessed in cross-validation inner loop
%		(default = 10)
%	'InnerPrintInfo': [true | false]
%		Print progress of cross-validation inner loop to the command window
%		(default = false)
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
%	'ErrBasedOn': ['unscaled', 'scaled']
%		Scaling of predictions, reconstructions and errors, whether they
%		should be returned and computed from scaled or unscaled entities (default
%		= 'unscaled'); note that PRESS, RMSEP, pred_bias and SEP will be
%		reported as scaled even in the unscaled errors are requested as they are
%		meant to assess the performance of the model and the comparison of errors
%		on different variables is possible only if they are all scaled
%	'ConLimCV': [lim]
%		Confidence limit for statistics on model parameters (default = 0.95)
%	'DOFMethodCV: ['naive' | 'cv_based' | 'krylov']
%		Method for computing degrees of freedom of the model, can be either based
%		on the number of latent variables ('naive', dof = A, default), based
%		on leave-one-out cross-validation (cv_based, dof estimated according
%		to Van Der Voet, 1999) or based on Krylov sequence ('krylov', dof
%		estimated according to Kramer, 2012)
%	'LimMethodCV': ['norm' | 't']
%		Method for computing the confidence limits of statsitcs on model
%		parameters, can be based on a normal distribution or on a t distribution
%		(default = 'norm')
%	'VariableSelection': ['none' | 'point' | 'pop']
%		Whether to perform predictor  variables selection according to VIP scores:
%		'none' (default) performs no predictor variable selection; 'point' keeps
%		only predictor variables such that the point estimate of the VIP is
%		greater than 1; 'pop' keeps only predictor variables such that the lower
%		confidence limit the VIP is greater than 1
%	'MultiYVarSel': ['any' | 'all']
%		Variable selection strategy for the multivariate Y case: 'any' (default)
%		keeps a predictor variable if it is significant for any of the response
%		variables; 'all' keeps a predictor variables only if it significant for
%		all of the response variables
%	'ConLimVIP': [lim]
%		Confidence limit for statistics on VIP scores for estimating their
%		confidence limits in cross-validation (default = 0.95)
%	'DOFMethodVIP': ['naive' | 'cv_based' | 'krylov']
%		Method for computing degrees of freedom of the model for estimation of
%		confidence limits on VIP scores and on estimates of control limits,
%		can be either based on the number of latent variables ('naive',
%		dof = A, default), based on leave-one-out cross-validation (cv_based,
%		dof estimated according to Van Der Voet, 1999) or based on Krylov
%		sequence ('krylov', dof estimated according to Kramer, 2012)
%	'LimMethodVIP': ['norm' | 't']
%		Method for computing the confidence limits of statistics on VIP scores
%		in cross-validation, can be based on a normal distribution or on a t
%		distribution (default = 'norm')
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
%	Note that variable seletion cannot be applied in a strightforward way to
%	nested cross validation, in particular different predictor variables could be
%	selected for different splits. However, estimates of error can still be be
%	obtained in order to comapre the prediction error of a model with selected
%	with variablee to tht error of a model with all variables. In this case,
%	prediction errors on different splits should be interpreted as different
%	estimates coming from hold-out validation splits, rather than from a
%	cross-validation loop.

%{
cross_validate_pls_nested.m
Version: 1.0.0
Date: 2021-10-06
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

cross_validate_pls_nested.m
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
method_inner = '';
G_obs_inner = [];
band_thick_obs_inner = [];
N_rep_inner = [];
A_inner = 10;
printinfo_inner = false;
preprocess = 'autoscale';
alg = 'nipals';
tol = 1e-11;
max_iter = 150;
norm_sch = 'standard';
err_on = 'unscaled';
varsel = 'none';
multysel = 'any';
predsel = false;
lim_VIP = 0.95;
dof_method_VIP = 'naive';
lim_method_VIP = 'norm';
lim_CV = 0.95;
dof_method_CV = 'naive';
lim_method_CV = 'norm';

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
			case 'InnerMethod'
				method_inner = varargin{i + 1};
			case 'InnerGroups'
				G_obs_inner = varargin{i + 1};
			case 'InnerBandThickness'
				band_thick_obs_inner = varargin{i + 1};
			case 'InnerRepetitions'
				N_rep_inner = varargin{i + 1};
			case 'InnerA'
				A_inner = varargin{i + 1};
			case 'InnerPrintInfo'
				printinfo_inner = varargin{i + 1};
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
			case 'VariableSelection'
				varsel = varargin{i + 1};
				if strcmp(varsel, 'none')
					predsel = false;
				elseif strcmp(varsel, 'point') || strcmp(varsel, 'pop')
					predsel = true;
				else
					error(['Supported variable selection methods are ' ...
						'none, point, and pop'])
				end
			case 'MultiYVarSel'
				multysel = varargin{i + 1};
				if strcmp(multysel, 'any') || strcmp(multysel, 'all')
				else
					error(['Supported variable selection methods for ' ...
						'multivariate Y are any and all'])
				end
			case 'ConLimVIP'
				lim_VIP = varargin{i + 1};
			case 'DOFMethodVIP'
				dof_method_VIP = varargin{i + 1};
				if strcmp(dof_method_VIP, 'naive') || strcmp(dof_method_VIP, 'cv_based') || strcmp(dof_method_VIP, 'krylov')
				else
					error(['Supported methods for degrees of freedom '...
						'estimation are naive, cv_based and krylov'])
				end
			case 'LimMethodVIP'
				lim_method_VIP = varargin{i + 1};
				if strcmp(lim_method_VIP, 'norm') || strcmp(lim_method_VIP, 't')
				else
					error(['Supported methods for VIP limits are normal '...
						'distribution method or t distribution method'])
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
				lim_method_CV = varargin{i + 1};
				if strcmp(lim_method_CV, 'norm') || strcmp(lim_method_CV, 't')
				else
					error(['Supported methods for limits are normal '...
						'distribution method or t distribution method'])
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

% Settings for inner loop
if isempty(method_inner)
	method_inner = 'continuous_blocks';
end
if isempty(G_obs_inner) && ~strcmp(method_inner, 'leave_one_out')
	if strcmp(method, 'leave_one_out')
		warning(['Better to specify the number of groups for the inner loop '...
			'when using leave_one_out in the outer loop: '...
			'setting G_obs_inner = 5 by default'])
		G_obs_inner = 5;
	else
		G_obs_inner = G_obs - 1;
	end
end
if isempty(band_thick_obs_inner)
	band_thick_obs_inner = 1;
end
if isempty(N_rep_inner)
	if ~strcmp(method_inner, 'random_subsets')
		N_rep_inner = 1;
	else
		error(['Number of inner repetitions must be provided using random '...
			'subsets as inner loop method'])
	end
end

% Discourage on using cv_based DOF for multivariate Y
if strcmp(dof_method_CV, 'cv_based') && V_Y > 1
	warning(['cv_based DOF estimation is meant for univariate Y: a variation'...
		'for multivariate Y will be used. However, for multivariate Y it is '...
		'suggested to use naive DOF estimation'])
end

%% Initialise results structure

NCV = struct;
	% Results
	NCV.results = struct;
		NCV.results.R_sq_Y = [];
		NCV.results.PRESS = [];
		NCV.results.RMSEP = [];
		NCV.results.pred_bias = [];
		NCV.results.SEP = [];
		NCV.results.A_sug = [];
		NCV.results.A_sug_freq = [];
		NCV.results.A_sug_pop = [];
	% Statistics on model parameters
	NCV.statistics = struct;
		NCV.statistics.R_sq_Y = [];
		NCV.statistics.PRESS = [];
		NCV.statistics.RMSEP = [];
		NCV.statistics.pred_bias = [];
		NCV.statistics.SEP = [];
		NCV.statistics.R = [];
	% Predictor variable selection
	NCV.selection = struct;
		NCV.selection.R_sq_Y_full = [];
		NCV.selection.PRESS_full = [];
		NCV.selection.RMSEP_full = [];
		NCV.selection.R_sq_Y_selected = [];
		NCV.selection.PRESS_selected = [];
		NCV.selection.RMSEP_selected = [];
		NCV.selection.VIP_mean = [];
		NCV.selection.VIP_cl = [];
		NCV.selection.selection = '';
		NCV.selection.multi_y = '';
		NCV.selection.lim_VIP = [];
		NCV.selection.dof_method_VIP = '';
		NCV.selection.lim_method_VIP = '';
	% Dimensions of the model entities
	NCV.dimensions = struct;
		NCV.dimensions.N = [];
		NCV.dimensions.V_X = [];
		NCV.dimensions.V_Y = [];
		NCV.dimensions.A = [];
		NCV.dimensions.N_rep = [];
	% Infos on cross-validation
	NCV.info = struct;
		NCV.info.CV_method_outer = '';
		NCV.info.grouping_outer = [];
		NCV.info.N_rep_outer = [];
		NCV.info.CV_method_inner = '';
		NCV.info.N_rep_inner = [];
		NCV.info.A_inner = [];
		NCV.info.lim = [];
		NCV.info.dof = [];
		NCV.info.pfac = [];
		NCV.info.preprocessing = '';
		NCV.info.normalisation = '';
		NCV.info.error_based_on = '';
		NCV.info.dof_method = '';
		NCV.info.obs_names = {};
		NCV.info.X_var_names = {};
		NCV.info.Y_var_names = {};

% Dimesions and labels assignments
NCV.dimensions.N = N;
NCV.dimensions.V_X = V_X;
NCV.dimensions.V_Y = V_Y;
NCV.dimensions.A = A;
NCV.dimensions.N_rep = N_rep;
NCV.info.lim = lim_CV;
NCV.info.pfac = pfac;
NCV.info.obs_names = obs_names;
NCV.info.X_var_names = X_var_names;
NCV.info.Y_var_names = Y_var_names;
NCV.info.preprocessing = preprocess;
NCV.info.normalisation = norm_sch;
NCV.info.error_based_on = err_on;
NCV.info.dof_method = dof_method_CV;

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
R_sq_pred_reps = zeros(A, V_Y, N_rep);
PRESS_reps = zeros(A, V_Y, N_rep);
pred_bias_reps = zeros(A, V_Y, N_rep);
SEP_reps = zeros(A, V_Y, N_rep);

if isempty(G_obs) && strcmp(method, 'leave_one_out')
	if aggr
		A_sug_reps = zeros(1, N, N_rep);
	else
		A_sug_reps = zeros(V_Y, N, N_rep);
	end
	R_sq_Y_split_full = zeros(A, V_Y, N, N_rep);
	PRESS_split_full = zeros(A, V_Y, N, N_rep);
	RMSEP_split_full = zeros(A, V_Y, N, N_rep);
	R_sq_Y_split_selected = zeros(A, V_Y, N, N_rep);
	PRESS_split_selected = zeros(A, V_Y, N, N_rep);
	RMSEP_split_selected = zeros(A, V_Y, N, N_rep);
	NCV.selection.VIP_mean = zeros(V_X, V_Y, A, N, N_rep);
	NCV.selection.VIP_cl = zeros(V_X, V_Y, A, N, N_rep);
else
	if aggr
		A_sug_reps = zeros(1, G_obs, N_rep);
	else
		A_sug_reps = zeros(V_Y, G_obs, N_rep);
	end
	R_sq_Y_split_full = zeros(A, V_Y, G_obs, N_rep);
	PRESS_split_full = zeros(A, V_Y, G_obs, N_rep);
	RMSEP_split_full = zeros(A, V_Y, G_obs, N_rep);
	R_sq_Y_split_selected = zeros(A, V_Y, G_obs, N_rep);
	PRESS_split_selected = zeros(A, V_Y, G_obs, N_rep);
	RMSEP_split_selected = zeros(A, V_Y, G_obs, N_rep);
	NCV.selection.VIP_mean = zeros(V_X, V_Y, A, G_obs, N_rep);
	NCV.selection.VIP_cl = zeros(V_X, V_Y, A, G_obs, N_rep);
end
R_pred_pop_reps = zeros(N, V_Y, A, N_rep);

grouping_reps = zeros(N, N_rep);

% Print info
if printinfo && strcmp(method, 'random_subsets')
	fprintf('Pre-allocations done, starting cross-validation outer loop\n')
	tic
end

% iteration loop
for i = 1:N_rep
	
	% Print info
	if printinfo && strcmp(method, 'random_subsets')
		fprintf('CV outer iteration: %.3g/%.3g\n', i, N_rep)
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
	R_pred = zeros(N, V_Y, A);
	R_pred_uns = zeros(N, V_Y, A);
	
	% Print info
	if printinfo && ~strcmp(method, 'random_subsets')
		fprintf('Pre-allocations done, starting cross-validation outer loop\n')
		tic
	end

	%% Cross-validation outer loop

	for g = 1:G_obs
		
		% Print progress of CV
		if printinfo && ~strcmp(method, 'random_subsets')
			fprintf('CV outer group: %.3g/%.3g\n', g, G_obs)
		end
		
		% Validation group index
		gti = g == grouping_obs;
		% Calibration dataset
		X_train = X(~gti, :);
		Y_train = Y(~gti, :);
		% Validation sample
		X_test = X(gti, :);
		Y_test = Y(gti, :);

		% Cross validation inner loop
		if strcmp(method_inner, 'random_subsets')
			switch alg
				case 'simpls'
					CV_inner = cross_validate_pls(X_train, Y_train, A_inner,...
						method_inner, G_obs_inner, N_rep_inner,...
						'Algorithm', alg,...
						'PrintInfo', printinfo_inner,...
						'Aggregate', aggr,...
						'ParFac', pfac,...
						'Preprocessing', preprocess,...
						'Normalise', norm_sch,...
						'VariableSelection', predsel,...
						'ConLimVIP', lim_VIP,...
						'DOFMethodVIP', dof_method_VIP,...
						'LimMethodVIP', lim_method_VIP...
					);
				case 'nipals'
					CV_inner = cross_validate_pls(X_train, Y_train, A_inner,...
						method_inner, G_obs_inner, N_rep_inner,...
						'Algorithm', alg, 'Tol', tol, 'MaxIter', max_iter,...
						'PrintInfo', printinfo_inner,...
						'Aggregate', aggr,...
						'ParFac', pfac,...
						'Preprocessing', preprocess,...
						'Normalise', norm_sch,...
						'VariableSelection', predsel,...
						'ConLimVIP', lim_VIP,...
						'DOFMethodVIP', dof_method_VIP,...
						'LimMethodVIP', lim_method_VIP....
					);
			end
		else
			switch alg
				case 'simpls'
					CV_inner = cross_validate_pls(X_train, Y_train, A_inner,...
						method_inner, G_obs_inner, band_thick_obs_inner,...
						'Algorithm', alg,...
						'PrintInfo', printinfo_inner,...
						'Aggregate', aggr,...
						'ParFac', pfac,...
						'Preprocessing', preprocess,...
						'Normalise', norm_sch,...
						'VariableSelection', predsel,...
						'ConLimVIP', lim_VIP,...
						'DOFMethodVIP', dof_method_VIP,...
						'LimMethodVIP', lim_method_VIP...
					);
				case 'nipals'
					CV_inner = cross_validate_pls(X_train, Y_train, A_inner,...
						method_inner, G_obs_inner, band_thick_obs_inner,...
						'Algorithm', alg, 'Tol', tol, 'MaxIter', max_iter,...
						'PrintInfo', printinfo_inner,...
						'Aggregate', aggr,...
						'ParFac', pfac,...
						'Preprocessing', preprocess,...
						'Normalise', norm_sch,...
						'VariableSelection', predsel,...
						'ConLimVIP', lim_VIP,...
						'DOFMethodVIP', dof_method_VIP,...
						'LimMethodVIP', lim_method_VIP...
					);
			end
		end
		
		%% Model building
		
		% Scale calibration data
		switch preprocess
			case 'none'
				mu_X = zeros(1, V_X);
				sigma_X = ones(1, V_X);
				X_train_s = X_train;
				mu_Y = zeros(1, V_Y);
				sigma_Y = ones(1, V_Y);
				Y_train_s = Y_train;
			case 'mean_centre'
				mu_X = mean(X_train);
				sigma_X = ones(1, V_X);
				X_train_s = X_train - mu_X;
				mu_Y = mean(Y_train);
				sigma_Y = ones(1, V_Y);
				Y_train_s = Y_train - mu_Y;
			case 'autoscale'
				[X_train_s, mu_X, sigma_X] = autoscale(X_train);
				[Y_train_s, mu_Y, sigma_Y] = autoscale(Y_train);
		end
		% Scale testing data
		X_test_s = scale_by(X_test, mu_X, sigma_X);
		Y_test_s = scale_by(Y_test, mu_Y, sigma_Y);
		
		% Build the model
		switch alg
			case 'simpls'
				[~, W_ast, ~, ~, Q, ~, b] = pls_calibration_simpls (X_train_s, Y_train_s, A, 'Normalise', norm_sch);
			case 'nipals'
				[~, W_ast, ~, ~, Q, ~, b] = pls_calibration_nipals(X_train_s, Y_train_s, A, 'Tol', tol, 'MaxIter', max_iter, 'Normalise', norm_sch);
		end
		
		% Compute scores applying the model
		T_test = X_test_s*W_ast;
		% Loop over LVs
		for a = 1:A
			% Residuals in cross-validation
			R_pred(gti, :, a) = Y_test_s - T_test(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
			% Unscaled cross-validation residuals for determination coefficients
			R_pred_uns(gti, :, a) = R_pred(gti, :, a)*diag(sigma_Y);
		end

		% Assignment to populations
		A_sug_reps(:, g, i) = CV_inner.results.A_sug';

		%% Variable selection

		% Check if variable selection is needed
		if predsel

			% Sum of squared Y for testing
			ssqY_test = sum((Y_test - mean(Y_test)).^2);
			% Number of observations in the current split
			N_g = sum(gti);
			% Loop over LVs
			for a = 1:A
				% Determination coefficient of full model
				R_sq_Y_split_full(a, :, g, i) = 1 - sum((R_pred_uns(gti, :, a)).^2, 1)./ssqY_test;
				% PRESS of full model
				PRESS_split_full(a, :, g, i) = sum((R_pred(gti, :, a)).^2, 1);
				% Root mean squared error in cross-validation
				RMSEP_split_full(a, :, g, i) = sqrt(PRESS_split_full(a, :, g, i)/N_g);
			end
		
			% Define variable selector
			switch varsel
				case 'none'
					% All predictors selected
					selector = CV_inner.selection.VIP_mean;
				case 'point'
					% Selector is mean of VIP scores
					selector = CV_inner.selection.VIP_mean;
				case 'pop'
					% Selector is lower confidence bound of VIP scores
					selector = CV_inner.selection.VIP_mean - CV_inner.selection.VIP_cl;
			end
			
			% Select variables
			switch multysel
				case 'any'
					% Variable is selected if significant for any response
					selected = any(permute(selector, [1, 3, 2]) > 1, 3);
				case 'all'
					% Variable is selected if significant for all responses
					selected = all(permute(selector, [1, 3, 2]) > 1, 3);
			end

			% Maxiumum number of LVs feasible
			A_feas = min(size(X_train_s(:, selected(:, a))));
			
			% Loop over feasible LVs
			for a = 1:A_feas
				% Check if any variable selected
				if any(selected(:, a))
					% Build the model
					switch alg
						case 'simpls'
							[~, W_ast, ~, ~, Q, ~, b] = pls_calibration_simpls (X_train_s(:, selected(:, a)), Y_train_s, a, 'Normalise', norm_sch);
						case 'nipals'
							[~, W_ast, ~, ~, Q, ~, b] = pls_calibration_nipals(X_train_s(:, selected(:, a)), Y_train_s, a, 'Tol', tol, 'MaxIter', max_iter, 'Normalise', norm_sch);
					end
			
					% Compute scores applying the model
					T_test = X_test_s(:, selected(:, a))*W_ast;
				
					% Residuals in cross-validation
					R_pred_selected = Y_test_s - T_test*diag(b)*Q';
					% Unscaled cross-validation residuals for determination coefficients
					R_pred_uns_selected = R_pred_selected*diag(sigma_Y);
					% Determination coefficient of full model
					R_sq_Y_split_selected(a, :, g, i) = 1 - sum((R_pred_uns_selected).^2, 1)./ssqY_test;
					% PRESS of full model
					PRESS_split_selected(a, :, g, i) = sum((R_pred_selected).^2, 1);
					% Root mean squared error in cross-validation
					RMSEP_split_selected(a, :, g, i) = sqrt(PRESS_split_selected(a, :, g, i)/N_g);
				else
					% No variable selected, assign NaNs
					R_sq_Y_split_selected(a, :, g, i) = NaN(1, V_Y);
					PRESS_split_selected(a, :, g, i) = NaN(1, V_Y);
					RMSEP_split_selected(a, :, g, i) = NaN(1, V_Y);
				end
			end
			% Loop over unfeasible LVs
			for a = (A_feas + 1):A
				% Infill other LVs with NaNs
				R_sq_Y_split_selected(a, :, g, i) = NaN(1, V_Y);
				PRESS_split_selected(a, :, g, i) = NaN(1, V_Y);
				RMSEP_split_selected(a, :, g, i) = NaN(1, V_Y);
			end
			
			% Assigments to the model structure
			NCV.selection.R_sq_Y_full = R_sq_Y_split_full;
			NCV.selection.PRESS_full = PRESS_split_full;
			NCV.selection.RMSEP_full = RMSEP_split_full;
			NCV.selection.R_sq_Y_selected = R_sq_Y_split_selected;
			NCV.selection.PRESS_selected = PRESS_split_selected;
			NCV.selection.RMSEP_selected = RMSEP_split_selected;
			NCV.selection.VIP_mean(:, :, :, g, i) =  CV_inner.selection.VIP_mean;
			NCV.selection.VIP_cl(:, :, :, g, i) = CV_inner.selection.VIP_cl;
			NCV.selection.selection = varsel;
			NCV.selection.multi_y = multysel;
			NCV.selection.lim_VIP = lim_VIP;
			NCV.selection.dof_method_VIP = dof_method_VIP;
			NCV.selection.lim_method_VIP = lim_method_VIP;
		else
			% Variable selection not requested
			NCV.selection.selection = varsel;
		end
	end
	
	% Determination coefficient in prediction
	R_sq_pred = 1 - permute(sum((R_pred_uns).^2, 1), [3, 2, 1])./sum((Y - mean(Y)).^2);
	% Prediction error sum of squares update
	PRESS = permute(sum((R_pred).^2, 1), [3, 2, 1]);
	% Prediction error bias
	pred_bias = permute(sum(R_pred, 1)/N, [3, 2, 1]);
	% Stadard error of prediction in prediction
	SEP = permute(sqrt(sum((R_pred - sum(R_pred, 1)/N).^2, 1)/(N - 1)), [3, 2, 1]);
	
	% Rescale entities according to preferences
	if strcmp(err_on, 'unscaled')
		R_pred = R_pred_uns;
	end
	
	% Save samples of populations of parameters in outer variables
	R_sq_pred_reps(:, :, i) = R_sq_pred;
	PRESS_reps(:, :, i) = PRESS;
	pred_bias_reps(:, :, i) = pred_bias;
	SEP_reps(:, :, i) = SEP;
	R_pred_pop_reps(:, :, :, i) = R_pred;
end

% Root mean squared error in cross-validation
RMSEP_reps = sqrt(PRESS_reps/N);

% Assigments to the model structure
NCV.info.CV_method_outer = method;
NCV.info.grouping_outer = grouping_reps;
NCV.info.N_rep_outer = N_rep;
NCV.info.dof = dof_CV;
NCV.info.CV_method_inner = method_inner;
NCV.info.N_rep_inner = N_rep_inner;
NCV.info.A_inner = A_inner;

%% Optimal number of LVs by inner loop

% Auxiliary index to avoid aggregation selector
auxi = size(A_sug_reps, 1);

% Find frequency of suggested number of LVs
A_sug_count = zeros(auxi, A);
% Loop over LVs
for a = 1:A
	% Count how many times a number of LVs has been suggested
	A_sug_count(:, a) = sum(reshape(A_sug_reps, auxi, G_obs*N_rep) == a, 2);
end
% Compute frquency
A_sug_freq = A_sug_count./sum(A_sug_count, 2);

% Final number of LVs (indexes of most frquently suggeste number of LVs)
[~, A_sug] = max(A_sug_freq, [], 2);

% Assigments to the model structure
NCV.results.A_sug = A_sug;
NCV.results.A_sug_pop = A_sug_reps;
NCV.results.A_sug_freq = A_sug_freq;

%% Aggregate over repetitions

% Note that no selection of the error to whole is explicitly done: performance
% indeces for LVs from 1 to A are returned in the model structure, together with
% the optimal values of A_sug. It is then upt to the user to explitly mark the
% performance for the number of LVs of interest.

% Determination coefficient in prediction
R_sq_pred = sum(R_sq_pred_reps, 3)/N_rep;
% Root mean squared error in prediction
PRESS = sum(PRESS_reps, 3)/N_rep;
% Root mean squared error in prediction
RMSEP = sum(RMSEP_reps, 3)/N_rep;
% Cross-validation error bias
pred_bias =  sum(pred_bias_reps, 3)/N_rep;
% Stadard error of prediction in prediction
SEP =  sum(SEP_reps, 3)/N_rep;

% Assigments to the model structure
NCV.results.R_sq_Y = R_sq_pred;
NCV.results.PRESS = PRESS;
NCV.results.RMSEP = RMSEP;
NCV.results.pred_bias = pred_bias;
NCV.results.SEP = SEP;

%% Confidence limits of performance indexes

% Print info
if printinfo
	fprintf('\nCross-validation outer loop done, estimating variability of parameters\n\n')
end

% Upper and lower confidence limits
% l_lim = (1 - lim_CV)/2;
u_lim = lim_CV + (1 - lim_CV)/2;

% Assigments to the model structure: poulations
NCV.statistics.R_sq_Y.pop = R_sq_pred_reps;
NCV.statistics.PRESS.pop = PRESS_reps;
NCV.statistics.RMSEP.pop = RMSEP_reps;
NCV.statistics.pred_bias.pop = pred_bias_reps;
NCV.statistics.SEP.pop = SEP_reps;
NCV.statistics.R.pop = R_pred_pop_reps;
NCV.statistics.R_sq_Y.mean = mean(R_sq_pred_reps, 3);
NCV.statistics.PRESS.mean = mean(PRESS_reps, 3);
NCV.statistics.RMSEP.mean = mean(RMSEP_reps, 3);
NCV.statistics.pred_bias.mean = mean(pred_bias_reps, 3);
NCV.statistics.SEP.mean = mean(SEP_reps, 3);
NCV.statistics.R.mean = mean(R_pred_pop_reps, 4);
% Assigments to the model structure: standard deviations
NCV.statistics.R_sq_Y.std = std(R_sq_pred_reps, 0, 3);
NCV.statistics.PRESS.std = std(PRESS_reps, 0, 3);
NCV.statistics.RMSEP.std = std(RMSEP_reps, 0, 3);
NCV.statistics.pred_bias.std = std(pred_bias_reps, 0, 3);
NCV.statistics.SEP.std = std(SEP_reps, 0, 3);
NCV.statistics.R.std = std(R_pred_pop_reps, 0, 4);
% Assigments to the model structure: limits around the mean
switch lim_method_CV
	case 'norm'
		% Confidence limits around a null value
		NCV.statistics.R_sq_Y.cl = norminv(u_lim, zeros(A, V_Y), NCV.statistics.R_sq_Y.std);
		NCV.statistics.PRESS.cl = norminv(u_lim, zeros(A, V_Y), NCV.statistics.PRESS.std);
		NCV.statistics.RMSEP.cl = norminv(u_lim, zeros(A, V_Y), NCV.statistics.RMSEP.std);
		NCV.statistics.pred_bias.cl = norminv(u_lim, zeros(A, V_Y), NCV.statistics.pred_bias.std);
		NCV.statistics.SEP.cl = norminv(u_lim, zeros(A, V_Y), NCV.statistics.SEP.std);
		NCV.statistics.R.cl = norminv(u_lim, zeros(N, V_Y, A), NCV.statistics.R.std);
	case 't'
		% Pre-allocations
		NCV.statistics.R_sq_Y.cl = zeros(A, V_Y);
		NCV.statistics.PRESS.cl = zeros(A, V_Y);
		NCV.statistics.RMSEP.cl = zeros(A, V_Y);
		NCV.statistics.pred_bias.cl = zeros(A, V_Y);
		NCV.statistics.SEP.cl = zeros(A, V_Y);
		NCV.statistics.R.cl = zeros(N, V_Y, A);
		% Loop over LVs
		for a = 1:A
			% Confidence limits around the mean
			DOF = max(N_rep - dof_CV(a), 0);
			t_cl = tinv(u_lim, DOF);
			NCV.statistics.R_sq_Y.cl(a, :) = sqrt(sum((NCV.statistics.R_sq_Y.pop(a, :, :) - NCV.statistics.R_sq_Y.mean(a, :)).^2, 3)/DOF)*t_cl;
			NCV.statistics.PRESS.cl(a, :) = sqrt(sum((NCV.statistics.PRESS.pop(a, :, :) - NCV.statistics.PRESS.mean(a, :)).^2, 3)/DOF)*t_cl;
			NCV.statistics.RMSEP.cl(a, :) = sqrt(sum((NCV.statistics.RMSEP.pop(a, :, :) - NCV.statistics.RMSEP.mean(a, :)).^2, 3)/DOF)*t_cl;
			NCV.statistics.pred_bias.cl(a, :) = sqrt(sum((NCV.statistics.pred_bias.pop(a, :, :) - NCV.statistics.pred_bias.mean(a, :)).^2, 3)/DOF)*t_cl;
			NCV.statistics.SEP.cl(a, :) = sqrt(sum((NCV.statistics.SEP.pop(a, :, :) - NCV.statistics.SEP.mean(a, :)).^2, 3)/DOF)*t_cl;
			NCV.statistics.R.cl(:, :, a) = sqrt(sum((NCV.statistics.R.pop(:, :, a, :) - NCV.statistics.R.mean(:, :, a)).^2, 4)/DOF)*t_cl;
		end
% 	case 'quantile'
% 		% Quantile-based confidence limits
% 		NCV.statistics.R_sq_Y.lcl = quantile(NCV.statistics.R_sq_Y.pop, l_lim, 3);
% 		NCV.statistics.PRESS.lcl = quantile(NCV.statistics.PRESS.pop, l_lim, 3);
% 		NCV.statistics.RMSEP.lcl = quantile(NCV.statistics.RMSEP.pop, l_lim, 3);
% 		NCV.statistics.pred_bias.lcl = quantile(NCV.statistics.pred_bias.pop, l_lim, 3);
% 		NCV.statistics.SEP.lcl = quantile(NCV.statistics.SEP.pop, l_lim, 3);
% 		NCV.statistics.R.lcl = quantile(NCV.statistics.R.pop, l_lim, 4);
% 		NCV.statistics.R_sq_Y.ucl = quantile(NCV.statistics.R_sq_Y.pop, u_lim, 3);
% 		NCV.statistics.PRESS.ucl = quantile(NCV.statistics.PRESS.pop, u_lim, 3);
% 		NCV.statistics.RMSEP.ucl = quantile(NCV.statistics.RMSEP.pop, u_lim, 3);
% 		NCV.statistics.pred_bias.ucl = quantile(NCV.statistics.pred_bias.pop, u_lim, 3);
% 		NCV.statistics.SEP.ucl = quantile(NCV.statistics.SEP.pop, u_lim, 3);
% 		NCV.statistics.R.ucl = quantile(NCV.statistics.R.pop, u_lim, 4);
end

% Print info
if printinfo
	t_ela = toc;
	fprintf('Cross-validation outer loop completed in %.4g second\n\n', t_ela)
end

%% Output assignments

NCV_out = NCV;

end