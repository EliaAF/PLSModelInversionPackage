function PRESS_out = slim_cv_pls (X_in, Y_in, A_in, varargin)
% SLIM_CV_PLS Slim version of PLS model cross-validation, computes just PRESS
%	using leave-one-out mode and in scaled error mode
%
% Syntax:
%	PRESS = slim_cv_pls (X, Y, A)
%	PRESS = slim_cv_pls (_, 'Key', value)
%
% Inputs:
%	X:			Array of predictor variables
%	Y:			Array of response variables
%	A:			Number of latent variables to assess in cross-validation
%
% Outputs:
%	PRESS:		Prediction error sum of squares
%
% Keys and Values:
%	'Preprocessing': ['none' | 'mean_centre' | 'autoscale']
%		Preprocessing method to be applied to the input data structures, can be:
%		no scaling at all ('none'); mean-centring only ('mean-centre');
%		mean-centring and scaling to unit variance ('autoscale', default); note
%		that the algorithm implemented by this function requires mean-centred
%		data matrices, therefore if the option 'none' is used the user is in
%		charge of providing mean-centred matrices, otherwise an error is issued
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

%{
slim_cv_pls.m
Version: 1.0.0
Date: 2021-08-28
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

slim_cv_pls.m
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
preprocess = 'autoscale';
alg = 'nipals';
tol = 1e-11;
max_iter = 150;
norm_sch = 'standard';

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
			otherwise
				error(['Key ' key ' undefined'])
		end
	end
end

%% Cross-validation loop

% Pre-allocations
R_CV = zeros(N, V_Y, A);

for n = 1:N
	% Calibration dataset
	X_cal = X;
	X_cal(n, :) = [];
	Y_cal = Y;
	Y_cal(n, :) = [];
	% Validation sample
	X_val = X(n, :);
	Y_val = Y(n, :);
	% Scale calibration data
	switch preprocess
			case 'none'
				mu_X = zeros(1, V_X);
				sigma_X = ones(1, V_X);
				X_cal_s = X_cal;
				mu_Y = zeros(1, V_Y);
				sigma_Y = ones(1, V_Y);
				Y_cal_s = Y_cal;
			case 'mean-centre'
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
			[~, W_ast, ~, ~, Q, ~, b] = pls_calibration_simpls (X_cal_s, Y_cal_s, A, 'Normalise', norm_sch);
		case 'nipals'
			[~, W_ast, ~, ~, Q, ~, b] = pls_calibration_nipals(X_cal_s, Y_cal_s, A, 'Tol', tol, 'MaxIter', max_iter, 'Normalise', norm_sch);
	end
	% Loop over LVs
	for a = 1:A
		% Apply the model
		Y_pred = X_val_s*W_ast(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
		% Residuals in cross-validation
		R_CV(n, :, a) = Y_val_s - Y_pred;
	end
end

% Prediction error sum of squares
PRESS = permute(sum((R_CV).^2, 1), [3, 2, 1]);

%% Output assignments

PRESS_out = PRESS;

end