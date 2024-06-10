function [B_out, W_ast_out, P_out, T_out, Q_out, U_out, b_out, W_out] = pls_calibration_nipals (X_in, Y_in, A_in, varargin)
% PLS_BY_NIPALS PLS model calibration using the non-linear iterative partial
%	least squares (NIPaLS) algorithm
%
% Syntax:
%	[B, W_ast, P, T, Q, U, b, W] = pls_calibration_nipals (X, Y, A)
%	[_] = pls_calibration_nipals (X, Y, A, 'Key', value)
%
% Inputs:
%	X:		Array of predictors variables to be used to calibrate the model
%	Y:		Array of response variables to be used to calibrate the model
%	A:		Number of latent variables
%
% Keys and Values:
%	'Tol': [tol]
%		Tolerance for convergence as 2-norm on relative scores variation
%		(default = 1e-15)
%	'MaxIter': [max_iter]
%		Maximum number of iterations allowed (default = 100)
%	'Normalise': ['standard' | 'loadings']
%		Normalisation scheme to be adopted, either no precise normalisation
%		('standard', default, as in the NIPaLS algorithm by Wold 2001) or both
%		the X loadings and Y loadings ('loadings', coherently with PCA as in
%		Geladi 1986)
%
% Outputs:
%	B:		Regression coefficient matrix
%	W_ast:	Corrected weights matrix
%	P:		X loadings matrix
%	T:		X scores matrix
%	Q:		Y loadings matrix
%	U:		Y scores matrix
%	b:		Inner-relation regression coefficients vector
%	W:		Weights matrix
%
%
% NOTE
%	The data arrays to be used to calibrate the model must be mean-centred (each
%	varible must have zero mean) and possibly scaled to unit variance (each
%	variable must have unit variance).
%
% NOTE
%	The key 'Normalise' governs the normalisation scheme applied to PLS model
%	entities. The value 'standard' is based on Wold 2001 and yields a model where
%	nothing is actually normalised and no inner-relation coefficients (b) are
%	provided; actually, they are provided as a unit vector. Exaplined variances
%	are not "stored" anywhere exactly. The 'loadings' options is based on Geladi
%	1986 and yields a model where X and Y loadings are normalised and the
%	variances are stored into X and Y scores respectively; note that these are
%	explained variances in reconstruction of X and Y data arrays, while the
%	explained variance on Y is here called predicted variance on Y (PV_Y) and it
%	is "stored" into the approximated Y score matrix obtained by means of the
%	inner-relation coefficients as T*diag(b).
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

%{
pls_calibration_nipals.m
Version: 1.0.0
Date: 2021-11-04
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_calibration_nipals.m
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

%% Optional arguments development

% Optionals initialization
tol = 1e-15;
max_iter = 100;
norm_sch = 'standard';
% Development cycle
if ~isempty(varargin)
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
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

%% Initial checks

% Check if the data array is unfolded
if size(X, 3) ~= 1 || size(Y, 3) ~= 1
	error('The data arrays must be unfolded for PLS model calibration')
end
% Check if the requested number of LVs is feasible
if A > min(size(X))
	error(['The number of latent variables cannot exceed min(size(X)) = '...
		num2str(min(size(X)))])
end
% Check if the data-array is mean-centered
if sum(mean(X) > 2e-8) ~= 0 || sum(mean(Y) > 2e-8) ~= 0
	error('The data arrays must be mean-centred for PLS model calibration')
end
% Check if there are missing values
if sum(ismissing(X), 'all') ~= 0 || sum(ismissing(Y), 'all') ~= 0
	error('Missing values found: PLS cannot be calibrated')
end

% Number of observations and number of variables
[N, V_X] = size(X);
V_Y = size(Y, 2);

%% PLS model building

% Initialisation
W = zeros(V_X, A);
P = zeros(V_X, A);
T = zeros(N, A);
Q = zeros(V_Y, A);
U = zeros(N, A);
b = zeros(1, A);
W_ast = zeros(V_X, A);
E = X;
R = Y;
G = eye(V_X, V_X);

% PLS by NIPaLS
switch norm_sch
	case 'standard'
		for a = 1:A
			u = ones(N, 1);
			t = u;
			iter = 0;
			err = tol + 1;
			while err >= tol && iter < max_iter
				iter = iter + 1;
				w = (E'*u)/(u'*u);
				w = w/norm(w);
				t_new = E*w; % = (E*w)/(w'*w); but w is normalised, denominator equals 1
				q = (R'*t_new)/(t_new'*t_new);
				u = (R*q)/(q'*q);
				err = norm(t_new - t)/norm(t_new);
				t = t_new;
			end
% 			fprintf('LV = %2d, Iter = %3d, err_abs = %.5e, err_rel = %.5e\n', a, iter, norm(t_new - t), err)
			if iter == max_iter
				warning(['Maximum number of iterations reached calculating LV ' num2str(a)])
			end
			p = (E'*t)/(t'*t);
			w_ast = G*w;
			T(:, a) = t;
			P(:, a) = p;
			U(:, a) = u;
			Q(:, a) = q;
			W(:, a) = w;
			W_ast(:, a) = w_ast;
			b(a) = 1;
			E = E - t*p';
			R = R - t*q';
			G = G - w_ast*p';
		end
	case 'loadings'
		for a = 1:A
			u = ones(N, 1);
			t = u;
			iter = 0;
			err = tol + 1;
			while err >= tol && iter < max_iter
				iter = iter + 1;
				w = (E'*u)/(u'*u);
				w = w/norm(w);
				t_new = E*w; % = (E*w)/(w'*w); but w is normalised, denominator equals 1
				q = (R'*t_new)/(t_new'*t_new);
				q = q/norm(q);
				u = R*q; % = (R*q)/(q'*q); but w is normalised, denominator equals 1
				err = norm(t_new - t)/norm(t_new);
				t = t_new;
			end
% 			fprintf('LV = %2d, Iter = %3d, err_abs = %.5e, err_rel = %.5e\n', a, iter, norm(t_new - t), err)
			if iter == max_iter
				warning(['Maximum number of iterations reached calculating LV ' num2str(a)])
			end
			p = (E'*t)/(t'*t);
			norm_p = norm(p);
			p = p/norm_p;
			t = t*norm_p;
			w = w*norm_p;
			w_ast = G*w;
			T(:, a) = t;
			P(:, a) = p;
			U(:, a) = u;
			Q(:, a) = q;
			W(:, a) = w;
			W_ast(:, a) = w_ast;
			b(a) = (u'*t)/(t'*t);
			E = E - t*p';
			R = R - t*b(a)*q';
			G = G - w_ast*p';
		end
end

% Regression coefficients
B = W_ast*diag(b)*Q';

% % Loadings directed as the largest component
% [~, X_index] = max(abs(P), [], 1);
% X_colsign = sign(P(X_index + (0:V_X:(A - 1)*V_X)));
% W = W.*X_colsign;
% P = P.*X_colsign;
% T = T.*X_colsign;
% W_ast = W_ast.*X_colsign;
% [~, Y_index] = max(abs(Q), [], 1);
% Y_colsign = sign(Q(Y_index + (0:V_Y:(A - 1)*V_Y)));
% Q = Q.*Y_colsign;
% U = U.*Y_colsign;
% b = b.*X_colsign.*Y_colsign;

%% Output assignments

B_out = B;
W_ast_out = W_ast;
P_out = P;
T_out = T;
Q_out = Q;
U_out = U;
b_out = b;
W_out = W;

end