function [dof_out, min_dof_out] = pls_dof(X_in, Y_in, A_in, varargin)
% PLS_DOF Degrees of freedom of PLS model
%
% Syntax:
%	dof = pls_dof (X, Y, A)
%	dof = pls_dof (X, Y, A, 'Key', value)
%	[dof, min_dof] = pls_dof (_)
%
% Inputs:
%	X:		Array of predictor variables
%	Y:		Array of response variables
%	A:		Number of latent variables over which compute dof
%
% Outputs:
%	dof:		Degrees of freedom of the PLS model form 1 to A LVs
%	min_dof:	Theoretical minimum degrees of freedom (only for Krylov method,
%				if computable, NaN othewise)
%
% Keys and Values:
%	'Method': ['naive' | 'cv_based' | 'krylov']
%		Method for computing degrees of freedom of the model, can be either based
%		on the number of latent variables ('naive', dof = A , default), based
%		on leave-one-out cross-validation (cv_based, dof estimated according
%		to Van Der Voet, 1999) or based on Krylov sequence ('krylov', dof
%		estimated according to Kramer, 2012)
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
%
%
% NOTE
%	Degrees of freedom are properly defined for PLS with univariate Y. This
%	functions handle the multivariate Y case by assuming that Y are independent,
%	hence returns a matrix of degrees of freedom containing values for each one
%	of the Y variables. If dof are estimate with the cv_based method, the
%	treatment is unclear as a single cross-validation cycle is used for error
%	estimation using multivariate Y models. In summary, degrees of freedom should
%	be trusted only in the univariate Y case, while care should be takes for the
%	multivariate case. Keeping this is mind, a possible option is to average over
%	dof of each Y variables, though no certainty of the reliability of the
%	results exists.
%
% NOTE
%	The estimation of degrees of freedom is reliable only if the preprocessing
%	method is set to 'autoscale'. More precisely, the variances of all variables
%	need to be as even as possible. Degrees of freedom will be computed also for
%	other pre-processing methods, but estimates may be questionable.
%
% NOTE
%	The algorithm for computing degrees of freedom based on Krylov sequence may
%	exhibit some numerical instability, returning negative values for a large
%	number of LVs. A threshold for the maximum number of degrees of freedom is
%	then automatically computed as the point where the ones computed using the
%	Krylov sequence equals the number of LVs plus one (minus a tolerance of 0.25)
%	for allowing some approximatione error). All the values after the critical
%	threshold are set to NaN.

%{
pls_dof.m
Version: 1.0.0
Date: 2021-09-05
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_dof.m
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
	error('Missing values found: estimated cannot be estimated')
end

% Number of observations and number of variables
N = size(X_unscaled, 1);
V_Y = size(Y_unscaled, 2);

%% Optional arguments development

% Optionals initialization
preprocess = 'autoscale';
alg = 'nipals';
tol = 1e-11;
max_iter = 150;
norm_sch = 'standard';
dof_method = 'naive';

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
			case 'Method'
				dof_method = varargin{i + 1};
				if strcmp(dof_method, 'naive') || strcmp(dof_method, 'cv_based') || strcmp(dof_method, 'krylov')
				else
					error(['Supported methods for degrees of freedom '...
						'estimation are naive, cv_based and krylov'])
				end
			otherwise
				error(['Key ' key ' undefined'])
		end
	end
end

%% Degrees of freedom

% Selection of the method for DOF estimation
switch dof_method
	case 'naive'
		% Naive DOF estimation based on number of LVs
		dof = repmat((1:A)', 1, V_Y);
		% No way of computing minimum DOF here
		min_dof = NaN;
	case 'cv_based'
		% Scale data
		switch preprocess
			case 'none'
				X = X_unscaled;
				Y = Y_unscaled;
			case 'mean-centre'
				X = X_unscaled - mean(X_unscaled);
				Y = Y_unscaled - mean(Y_unscaled);
			case 'autoscale'
				X = autoscale(X_unscaled);
				Y = autoscale(Y_unscaled);
		end
		
		% Build the model
		switch alg
			case 'simpls'
				[~, ~, ~, T, Q, ~, b] = pls_calibration_simpls (X, Y, A, 'Normalise', norm_sch);
			case 'nipals'
				[~, ~, ~, T, Q, ~, b, ~] = pls_calibration_nipals(X, Y, A, 'Tol', tol, 'MaxIter', max_iter, 'Normalise', norm_sch);
		end
		
		% Sum of squared errors for all LVs
		SSE = zeros(A, V_Y);
		for a = 1:A
			% Prediction error with LV up to a
			R_a = Y - T(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
			% SSE in calibration (prediction error on Y)
			SSE(a, :) = sum(R_a.^2, 1);
		end
		% MSE in prediction based on redistribution (MSE in calibration)
		MSEP_RS = SSE/N;
		
		% Prediction error sum of squares from leave-one-out CV (slimmed)
		PRESS = slim_cv_pls(X_unscaled, Y_unscaled, A,...
			'Preprocessing', preprocess,...
			'Algorithm', alg,...
			'Tol', tol,...
			'MaxIter', max_iter,...
			'Normalise', norm_sch);
		
		% MSE in prediction based on CV
		MSEP_CV = PRESS/N;
		
		% Pseudo-degrees of freedom of PLS model
		dof = N*(1 - 1/N - sqrt(1/N^2 + (MSEP_RS./MSEP_CV)*(1 - 2/N - 2/(N*(N - 1)))));
% 		% The more general equation is being used, the simplified one is below
% 		dof = N*(1 - sqrt(MSEP_RS./MSEP_CV));

		% No way of computing minimum DOF here
		min_dof = NaN;
	case 'krylov'
		% Choiche of the preprocessing method
		switch preprocess
			case 'none'
				X = X_unscaled;
				Y = Y_unscaled;
			case 'mean_centre'
				X = X_unscaled - mean(X_unscaled);
				Y = Y_unscaled - mean(Y_unscaled);
			case 'autoscale'
				X = autoscale(X_unscaled);
				Y = autoscale(Y_unscaled);
		end
		
		% Selection of the algorithm for model calibration
		switch alg
			case 'simpls'
				[~, ~, ~, T, Q, ~, b] = pls_calibration_simpls (X, Y, A, 'Normalise', norm_sch);
			case 'nipals'
				[~, ~, ~, T, Q, ~, b, ~] = pls_calibration_nipals(X, Y, A, 'Tol', tol, 'MaxIter', max_iter, 'Normalise', norm_sch);
		end
		
		% Residuals for all LVs
		Y_pred = zeros(N, V_Y, A);
		switch preprocess
			case 'none'
				for a = 1:A
					Y_pred(:, :, a) = T(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
				end
			case 'mean_centre'
				sigma_Y = var(Y_unscaled);
				for a = 1:A
					Y_pred(:, :, a) = rescale_by(T(:, 1:a)*diag(b(1:a))*Q(:, 1:a)', zeros(1, V_Y), sigma_Y);
				end
			case 'autoscale'
				for a = 1:A
					Y_pred(:, :, a) = T(:, 1:a)*diag(b(1:a))*Q(:, 1:a)';
				end	
		end
		
		% Normalise T
		T = T/diag(sqrt(sum(T.^2, 1)));

		% Covariance of X
		S = X'*X/(N - 1);
		% Maximum eigenvalue of covariance of X
		lambda = eigs(S, 1);
		if lambda <= 0.5*trace(S)
			% Minimum DOF
			min_dof = 1 + trace(S)/lambda;
		else
			% No way of computing minimum DOF
			min_dof = NaN;
		end
		
		% Kernel matrix of X
		K = X*X';

		% Eigenvalues of the kernel matrix in descending order
		lambda = sort(eig(K), 'descend');
		% Remove insignficant eigenvalues (tolerance adapted from pinv)
		lambda(lambda < N*eps(norm(lambda, inf))) = [];
		% Initialise vector for trace of the kernel matrix raised to a
		trK = zeros(A, 1);
		% Loop over latent variables
		for a = 1:A
			% Trace of the kernel matrix raised to a
			trK(a) = sum(lambda.^a);
		end

		% Initialise tensore for computing t_a*K^a*T as t_a*KaT
		KaT = zeros(N, A, A);
		% First element of KaT
		KaT(:, :, 1) = K*T;
		% Other elements by recursion
		for a = 2:A
			KaT(:, :, a) = K*KaT(:, :, a - 1);
		end

		% Auxiliary variable for computing c as inv(B)*b = inv(B)*T'*y
		b = T'*Y;

		% Pre-allocations
		dof = zeros(A, V_Y);
		% Loop over quality variables
		for w = 1:V_Y
			% Krylov sequence of K and y
			KY = zeros(N, A);
			% First element of the sequence
			KY(:, 1) = K*Y(:, w);
			% Elements of the sequence by recursion
			for a = 2:A
				KY(:, a) = K*KY(:, a - 1);
			end

			% Basis transformation to improve numerical stability
			B = T'*KY;
			% Trash useless components in B to improve numerical stability
			B = triu(B);
			% Invert the basis transformation
			warning('off', 'MATLAB:nearlySingularMatrix')
			B_inv = B\eye(A);
			warning('on', 'MATLAB:nearlySingularMatrix')
			% Note: B is singular and MatLab will issues a warning. However, if you
			% regularise the inversion using piv, dof falls on the naive ones, interestingly.
			% to still allow for the inversion of the singular matrix yet avoiding the
			% warning, disable it before inverting and re-enable it afterwards.

			% Pre-allocation of terms for dof
			ctrK_a = zeros(A, 1);	% Trace term
			tKat = zeros(A, 1);		% Kernel term
			yKav = zeros(A, 1);		% Residual term

			% Loop over latent variables
			for a = 1:A
				% Entites up to a LVs
				B_inv_a = B_inv(1:a, 1:a);
				c_a = B_inv_a*b(1:a, w);
				V_a = T(:, 1:a)*B_inv_a';
				R_a = Y(:, w) - Y_pred(:, w, a);
				% Trace term
				ctrK_a(a) = c_a'*trK(1:a);
				% Loop over latent variables
				for j = 1:a
					% Kernel term
					tKat(a) = tKat(a) + c_a(j)*trace(T(:, 1:a)'*KaT(:, 1:a, j));
					% Recursion of residuals
					R_a = K*R_a;
					% Residual term
					yKav(a) = yKav(a) + R_a'*V_a(:, j);
				end
			end

			% Degrees of freedom: trace term - kernel term + residual term + LVs + intercept
			dof(:, w) = ctrK_a - tKat + yKav + (1:A)' + 1;
			
			% Check for maximum number of DOFs
			idx = dof(:, w) < (1:A)' + 1 - 0.25;
			% Set to NaN meaningless values
			dof(idx, w) = NaN;
		end
end

%% Output assignments

dof_out = dof;
min_dof_out = min_dof;

end