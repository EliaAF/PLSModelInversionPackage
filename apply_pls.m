function pred_out = apply_pls (model_in, X_in, Y_in, varargin)
% APPLY_PLS PLS model application for prediction
%
% Syntax:
%	pred = apply_pls (model, X, Y)
%	pred = apply_pls (model, X, Y, 'Key', value)
%
% Inputs:
%	model:	Calibrated PLS model (see build_pls function)
%	X:		Array of predictors variables to apply the model to
%	Y:		Array of response variables to be used for error estimation (if not
%			available set it to [], corresponding errors and diagnostics will be
%			populated by NaNs)
%
% Outputs:
%	pred:	Structure of the prediction
%
% Keys and Values:
%	'ErrBasedOn': ['unscaled', 'scaled']
%		Scaling of predictions, reconstructions and errors, whether they
%		should be returned and computed from scaled or unscaled entities
%		Kind of contributions to diagnostics, which can be simple indictors
%		either positive or negative according to the direction of deviation and
%		with two-tail confidence limits (approach of Miller, 1998), or purely
%		positive contributions (like "squadred") that sum up the diagnostic
%		values and with one-tail confidence limits (approach of Westerhuis, 2000)
%		(default = 'simple')
%	'ConLim': [lim]
%		Confidence limit
%	'DOFMethod: ['naive' | 'cv_based' | 'krylov']
%		Method for computing degrees of freedom of the model, can be either based
%		on the number of latent variables ('naive', dof = A, default), based
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
%		plot, compute from a F distribution with considering only two principal
%		components (A_F = 2) or all the requested ones (A_F = A, default)
%	'ObsNames': [obs_names]
%		Names of the observations as chars in a cell array (default are
%		progressive numerical identifiers prefixed by the letter O)
%
%
% NOTE: Default values of the optinal arguments are inherithed from the PLS model
%	passed as input

%{
apply_pls.m
Version: 1.0.1
Date: 2022-12-20
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

apply_pls.m
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

model = model_in;
X_unscaled = X_in;
Y_unscaled = Y_in;

%% Initial checks

% Check if there are missing values
if sum(ismissing(X_unscaled), 'all') ~= 0 || sum(ismissing(Y_unscaled), 'all') ~= 0
	error('Missing values found: PLS cannot be applied')
end
% Check if the Y matrix has been provided
if isempty(Y_unscaled)
	% Make a fake Y matrix filled with NaNs
	Y_unscaled = NaN(size(X_unscaled, 1), model.dimensions.V_Y);
end
% Check if the data array is unfolded
if size(X_unscaled, 3) ~= 1 || size(Y_unscaled, 3) ~= 1
	error('The data arrays must be unfolded for PLS model calibration')
end
% Check if X variables match
if size(X_unscaled, 2) ~= model.dimensions.V_X
	error('Incorrect number of X variables')
end
% Check if Y variables match
if size(Y_unscaled, 2) ~= model.dimensions.V_Y
	error('Incorrect number of Y variables')
end

% Number of observations and number of variables
[N, V_X] = size(X_unscaled);
V_Y = size(Y_unscaled, 2);

%% Optional arguments development

% Optionals initialization
err_on = model.info.error_based_on;
lim = model.estimates.lim;
dof_method = model.info.dof_method;
Tsq_lim_method = model.info.Tsq_lim_method;
SqE_lim_method = model.info.SqE_lim_method;
con_lim_method = model.info.con_lim_method;
l_kind = model.info.l_kind;

obs_names = cellstr([repmat('O', N, 1) char(pad(replace(string(num2str((1:N)')),  ' ', ''), length(num2str(N)), 'left', '0')) repmat('_val', N, 1)]);

% Development cycle
if ~isempty(varargin)
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
			case 'ErrBasedOn'
				err_on = varargin{i + 1};
				if strcmp(err_on, 'scaled') || strcmp(err_on, 'unscaled')
				else
					error('Undefined key for ErrBasedOn')
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
					error(['Supported methods for limit on T_sq are '...
						'chi^2 method, F distribution method and '...
						'kernel density estimation method'])
				end
			case 'SqErrLimMethod'
				SqE_lim_method = varargin{i + 1};
				if strcmp(SqE_lim_method, 'chisq') || strcmp(SqE_lim_method, 'jack_mud') || strcmp(SqE_lim_method, 'KDE')
				else
					error(['Supported methods for limit on SRE_X are '...
						'chi^2 method, Jackson-Mudholkar equation and '...
						'kernel density estimation method'])
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
			otherwise
				error(['Key ' key ' undefined'])
		end
	end
end

%% Model structure initialization

% Model structure
pred = struct;
	% Scores, predictions and errors
	pred.prediction = struct;
		pred.prediction.T = [];
		pred.prediction.U = [];
		pred.prediction.Y_pred = [];
		pred.prediction.X_rec = [];
		pred.prediction.Y_rec = [];
		pred.prediction.R = [];
		pred.prediction.E = [];
		pred.prediction.F = [];
	% Data used for validation
	pred.data = struct;
		pred.data.X = [];
		pred.data.Y = [];
		pred.data.X_uns = [];
		pred.data.Y_uns = [];
	% Dimensions of the model entities
	pred.dimensions = struct;
		pred.dimensions.N = [];
		pred.dimensions.V_X = [];
		pred.dimensions.V_Y = [];
		pred.dimensions.A = [];
	% Sclaing applied to data
	pred.scaling = struct;
		pred.scaling.mu_X = [];
		pred.scaling.sigma_X = [];
		pred.scaling.mu_Y = [];
		pred.scaling.sigma_Y = [];
	% Perfomance of the model
	pred.performance = struct;
		pred.performance.R_sq_Y = [];
		pred.performance.SSE = [];
		pred.performance.RMSE = [];
		pred.performance.bias = [];
		pred.performance.SEP = [];
	% Diagnostics on the model
	pred.diagnostics = struct;
		pred.diagnostics.T_sq_T = [];
		pred.diagnostics.T_sq_U = [];
		pred.diagnostics.SPE_Y = [];
		pred.diagnostics.SRE_X = [];
		pred.diagnostics.SRE_Y = [];
		pred.diagnostics.T_sq_T_con = [];
		pred.diagnostics.T_sq_U_con = [];
		pred.diagnostics.T_sq_T_con_dir = [];
		pred.diagnostics.SPE_Y_con = [];
		pred.diagnostics.SRE_X_con = [];
		pred.diagnostics.SRE_Y_con = [];
	% Confidence limits
	pred.estimates = struct;
		pred.estimates.lim = [];
		pred.estimates.dof = [];
		pred.estimates.lim_T_sq_T = [];
		pred.estimates.lim_T_sq_U = [];
		pred.estimates.lim_SPE_Y = [];
		pred.estimates.lim_SRE_X = [];
		pred.estimates.lim_SRE_Y = [];
		pred.estimates.lim_T_sq_T_con = [];
		pred.estimates.lim_T_sq_U_con = [];
		pred.estimates.lim_T_sq_T_con_dir = [];
		pred.estimates.lim_SPE_Y_con = [];
		pred.estimates.lim_SRE_X_con = [];
		pred.estimates.lim_SRE_Y_con = [];
		pred.estimates.l_T = [];
		pred.estimates.l_U = [];
	% Infos on the model
	pred.info = struct;
		pred.info.obs_names = {};
		pred.info.X_var_names = {};
		pred.info.Y_var_names = {};
		pred.info.preprocessing = '';
		pred.info.error_based_on = '';
		pred.info.contribution_method = '';
		pred.info.dof_method = '';
		pred.info.Tsq_lim_method = '';
		pred.info.SqE_lim_method = '';
		pred.info.con_lim_method = '';
		pred.info.l_kind = '';

% Dimesions and labels assignments
pred.dimensions.N = N;
pred.dimensions.V_X = V_X;
pred.dimensions.V_Y = V_Y;
pred.info.obs_names = obs_names;
pred.info.X_var_names = model.info.X_var_names;
pred.info.Y_var_names = model.info.Y_var_names;
pred.info.preprocessing = model.info.preprocessing;
pred.info.error_based_on = err_on;
pred.info.contribution_method = model.info.contribution_method;
pred.info.dof_method = dof_method;
pred.info.Tsq_lim_method = Tsq_lim_method;
pred.info.SqE_lim_method = SqE_lim_method;
pred.info.con_lim_method = con_lim_method;
pred.info.l_kind = l_kind;

%% Entities from the model

A = model.dimensions.A;

mu_X = model.scaling.mu_X;
sigma_X = model.scaling.sigma_X;
mu_Y = model.scaling.mu_Y;
sigma_Y = model.scaling.sigma_Y;

B = model.parameters.B;
W_ast = model.parameters.W_ast;
Q = model.parameters.Q;
P = model.parameters.P;

sigma_sq_T = model.parameters.sigma_sq_T;
sigma_sq_U = model.parameters.sigma_sq_U;

% Assigments to the model structure
pred.dimensions.A = A;

%% Preprocessing

% Apply the preprocessing method
X = scale_by(X_unscaled, mu_X, sigma_X);
Y = scale_by(Y_unscaled, mu_Y, sigma_Y);

% Assigments to the model structure
pred.data.X = X;
pred.data.Y = Y;
pred.data.X_uns = X_unscaled;
pred.data.Y_uns = Y_unscaled;
pred.scaling.mu_X = mu_X;
pred.scaling.sigma_X = sigma_X;
pred.scaling.mu_Y = mu_Y;
pred.scaling.sigma_Y = sigma_Y;

%% Model application and prediction performance

% Predicted Y
Y_pred = X*B;
% Predicted scores
T = X*W_ast;
U = Y*Q*pinv(Q'*Q);
% Reconstructued matrices
X_rec = T*P';
Y_rec = U*Q';
% Y prediction residuals
R = Y - Y_pred;
% Reconstruction errors
E = X - X_rec;
F = Y - Y_rec;

% Determination coefficient on each response variable
R_sq_Y = 1 - diag(R'*R)./diag(Y'*Y);
% SSE in prediction (prediction error on Y)
SSE = sum(R.^2, 1);
% Root mean squared error
RMSE = sqrt(SSE/N);
% Prediction error bias
bias = sum(R, 1)/N;
% Standard error of prediction
SEP = sqrt(sum((R - bias).^2, 1)/(N - 1));

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
pred.prediction.T = T;
pred.prediction.U = U;
pred.prediction.Y_pred = Y_pred;
pred.prediction.X_rec = X_rec;
pred.prediction.Y_rec = Y_rec;
pred.prediction.R = R;
pred.prediction.E = E;
pred.prediction.F = F;
pred.performance.R_sq_Y = R_sq_Y;
pred.performance.SSE = SSE;
pred.performance.RMSE = RMSE;
pred.performance.bias = bias;
pred.performance.SEP = SEP;

%% Prediction diagnostics

% Hotelling's T^2 of T and U
T_sq_T = T_sq_statistic(T, sigma_sq_T);
T_sq_U = T_sq_statistic(U, sigma_sq_U);
% Squared prediction error (SPE) on Y
SPE_Y = Q_statistic(R_fd);
% Squared reconstruction errors (SRE) of X and Y
SRE_X = Q_statistic(E_fd);
SRE_Y = Q_statistic(F_fd);

% Choice of the kind of contributions
switch model.info.contribution_method
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
		% Contributions to SRE onA X and Y
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

% Assigments to the model structure
pred.diagnostics.T_sq_T = T_sq_T;
pred.diagnostics.T_sq_U = T_sq_U;
pred.diagnostics.SPE_Y = SPE_Y;
pred.diagnostics.SRE_X = SRE_X;
pred.diagnostics.SRE_Y = SRE_Y;
pred.diagnostics.T_sq_T_con = T_sq_T_con;
pred.diagnostics.T_sq_U_con = T_sq_U_con;
pred.diagnostics.T_sq_T_con_dir = T_sq_T_con_dir;
pred.diagnostics.SPE_Y_con = SPE_Y_con;
pred.diagnostics.SRE_X_con = SRE_X_con;
pred.diagnostics.SRE_Y_con = SRE_Y_con;

%% Estimation of confidence limits on validation

% Degrees of freedom of the model
dof = model.estimates.dof;

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
			warning(['Confidence ellipse requested as based on two '...
				'components only, possible only using F distribution-based '...
				'confidence limits for T_sq_T and T_sq_U: limits will be '...
				'inconsistent'])
		end
		lim_T_sq_T_for_plots = (2*(N - 1)/(N - 2))*finv(lim, 2, N - 2);
		lim_T_sq_U_for_plots = lim_T_sq_T_for_plots;
	case 'full'
		lim_T_sq_T_for_plots = lim_T_sq_T;
		lim_T_sq_U_for_plots = lim_T_sq_U;
end
% Confidence ellipsoid semiaxes of T and U
l_T = sqrt(sigma_sq_T*lim_T_sq_T_for_plots);
l_U = sqrt(sigma_sq_U*lim_T_sq_U_for_plots);

% Assigments to the model structure
pred.estimates.lim = lim;
pred.estimates.dof = dof;
pred.estimates.lim_T_sq_T = lim_T_sq_T;
pred.estimates.lim_T_sq_U = lim_T_sq_U;
pred.estimates.lim_SPE_Y = lim_SPE_Y;
pred.estimates.lim_SRE_X = lim_SRE_X;
pred.estimates.lim_SRE_Y = lim_SRE_Y;
pred.estimates.lim_T_sq_T_con = lim_T_sq_T_con;
pred.estimates.lim_T_sq_U_con = lim_T_sq_U_con;
pred.estimates.lim_T_sq_T_con_dir = lim_T_sq_T_con_dir;
pred.estimates.lim_SPE_Y_con = lim_SPE_Y_con;
pred.estimates.lim_SRE_X_con = lim_SRE_X_con;
pred.estimates.lim_SRE_Y_con = lim_SRE_Y_con;
pred.estimates.l_T = l_T;
pred.estimates.l_U = l_U;

%% Output assignments

pred_out = pred;

end