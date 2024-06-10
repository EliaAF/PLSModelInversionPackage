function Y_pred_cl_out = pls_prediction_uncertainty (T_pred_in, alpha_in, sigma_sq_in, SSE_in, dof_in, N_in, varargin)
% PLS_PREDICTION_UNCERTAINTY Uncertainty on prediction of PLS model
%
% Syntax:
%	Y_pred_cl = pls_prediction_uncertainty (T_pred, alpha, sigma_sq, SSE, dof, N)
%	Y_pred_cl = pls_prediction_uncertainty (T_pred, alpha, sigma_sq, SSE, dof, N, 'Key', value)
%
% Inputs:
%	T_pred:		Scores matrix of the predicted values for uncertainty estimation
%	alpha:		Significance level of the confidence limit
%	sigma_sq:	Variances of the latent variables in the calibration dataset
%	SSE:		Sum of squared errors of the output variables in calibration
%	dof:		Degrees of freedom of the PLS model
%	N:			Number of observations in the calibration dataset of the PLS
%				model
%
% Outputs:
%	Y_pred_cl:	Confidence limit of the predictions
%
% Keys and Values:
%	'Limit': ['prediction', 'confidence']
%		Whther the uncertainty should be evaluated according to the prediction
%		limit approach (default) or based on the confidence deriving from the
%		calibration dataset
%
%
% NOTE
%	The SSE must be provided as a row vector contianing SSEs for the single
%	output variables of the PLS model (estimated from the calibration data). In
%	the same way, the confidence limit is provided as a matrix reporting the
%	limit for every observation in T_new for all output variables. For example,
%	if T_new has been used to estimate Y_pred, then the uncertinaty is defined
%	as:
%		Y_pred +- Y_pred_cl
%	The confidence limit is returned according the the two tail test, its value
%	scaled on the standard deviations of the output vairables in the calibration
%	dataset.
%
% NOTE
%	By default, the function yields a prediction limit. This limit accounts for
%	an "additional" standard error to guarantee that the confidence region around
%	the point estimate (the predicted value) of the output variable covers the
%	true value at significance level alpha. On the other hand, the simple
%	confidence limit offers this guarantee only around the expected value
%	(average of the distribution) of the model prediction, thus being less
%	conservative.

%{
pls_prediction_uncertainty.m
Version: 1.0.0
Date: 2023-12-21
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_prediction_uncertainty.m
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

T_pred = T_pred_in;
alpha = alpha_in;
sigma_sq = sigma_sq_in;
SSE = SSE_in;
dof = dof_in;
N = N_in;

%% Optional arguments development

% Optionals initialization
test_cl = 'prediction';

% Development cycle
if ~isempty(varargin)
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
			case 'Limit'
				test_cl = varargin{i + 1};
				if strcmp(test_cl, 'prediction') || strcmp(test_cl, 'confidence')
				else
					error('Supported limits are prediction and confidence')
				end
			otherwise
				error(['Key ' key ' undefined'])
		end
	end
end

%% PLS model prediction uncertainty

% Determine limit factor
switch test_cl
	case 'prediction'
		% Compute prediction limit
		lf = 1;
	case 'confidence'
		% Compute simple confidence limit
		lf = 0;
end

% Limit from significance level
lim = 1 - alpha/2;

% Root-mean-squared error of the calibration dataset (DOF-corrected)
RMSE = sqrt(SSE/(N - dof));
% Leverages of predictions
h_new = sum((T_pred.^2)*diag(sigma_sq.^-1), 2)/(N - 1);
% Standard deviation of prediction
s_Y_pred = RMSE.*sqrt(lf + 1/N + h_new);

% Confidence limit for the prediction
Y_pred_cl = s_Y_pred*tinv(lim, N - dof);

%% Output assignments

Y_pred_cl_out = Y_pred_cl;

end