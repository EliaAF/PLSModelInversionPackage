function [B_out, W_ast_out, P_out, T_out, Q_out, U_out, b_out, V_out] = pls_calibration_simpls (X_in, Y_in, A_in, varargin)
% PLS_BY_SIMPLS PLS model calibration using the statistically-inspired
%	modification to partial least squares (SIMPLS) algorithm
%
% Syntax:
%	[B, W_ast, P, T, Q, U, b, V] = pls_calibration_simpls (X, Y, A)
%	[_] = pls_calibration_simpls (X, Y, A, 'Key', value)
%
% Inputs:
%	X:		Array of predictor variables to be used to calibrate the model
%	Y:		Array of response variabls to be used to calibrate the model
%	A:		Number of latent variables
%
% Outputs:
%	B:		Regression coefficient matrix
%	W_ast:	Corrected weights matrix
%	P:		X loadings matrix
%	T:		X scores matrix
%	Q:		Y loadings matrix
%	U:		Y scores matrix
%	b:		Inner-relation regression coefficients vector
%	V:		Basis of the X loadings space
%
% Keys and Values:
%	'Normalise': ['standard' | 'loadings']
%		Normalisation scheme to be adopted, either to normalise X scores
%		('standard', default, as in the SIMPLS algorithm by de Jong 1993) or both
%		the X loadings and Y loadings ('loadings', coherently with PCA as in
%		Geladi 1986)
%
%
% NOTE
%	The data arrays to be used to calibrate the model must be mean-centred (each
%	varible must have zero mean) and possibly scaled to unit variance (each
%	variable must have unit variance).
%
% NOTE
%	The key 'Normalise' governs the normalisation scheme applied to PLS model
%	entities. The value 'standard' is based on de Jong 1993 and yields a model
%	where the X scores are normalised and no inner-relation coefficients (b) are
%	provided. Variances of X and Y are stored into X and Y loadings respectively;
%	note that the variance on X is related to the reconstruction of X, therefore
%	it is called explained variance on X (EV_X), while variance on Y is related
%	to the prediction of Y, therefore it is called predicted variance on Y
%	(PV_Y). There is no way for getting the explained variance in reconstruction
%	of Y without modifying the normalisation scheme of the original algorithm. In
%	fact, this function implements a slight variation of the normalisations
%	scheme where the Y scores are normalised for this reason. the 'loadings'
%	options is based on the scaling proposed by Geladi 1986 for the NIPaLS
%	algorithm and yields a model where X and Y loadings are normalised and the
%	variances are stored in X and Y scores respectively; as for the 'standard'
%	value, these variances are EV_X and PV_Y respectively; one can also get EV_Y
%	with this normalisation scheme, but this is just an approximation of the
%	"true" reconstruction variance of Y: see the papaer by de Jong for details.
%	Finally, note that inner-relation coefficients are not provided in the
%	original paper, while here they are provided as a unit vector with the
%	'standard' option and as real coefficients in with the 'loadings' option.
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
%	The workflow for prediction is a three-step procedure:
%		1) Decompose a new sample x using the corrected weights:
%			t = x*W_ast
%		2) Regress the Y scores from the X scores by means of the inner-relation
%		   regression coefficients:
%			u_approx = t*diag(b)
%		3) Compose the prediction of the y using the Y loadings:
%			y_approx = u_approx*Q
%	The workflow can be summarised into a single operation defining the PLS
%	regression coefficients matrix:
%		B = W_ast*diag(b)*Q'
%	thus:
%		y_approx = x*B

%{
pls_calibration_simpls.m
Version: 1.0.0
Date: 2021-09-24
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_calibration_simpls.m
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
norm_sch = 'standard';
% Development cycle
if ~isempty(varargin)
	for i = 1:2:length(varargin)
		key = varargin{i};
		switch key
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
if sum(mean(X) > 1e-9) ~= 0 || sum(mean(Y) > 1e-9) ~= 0
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
W_ast = zeros(V_X, A);
P = zeros(V_X, A);
T = zeros(N, A);
Q = zeros(V_Y, A);
U = zeros(N, A);
V = zeros(V_X, A);
S = X'*Y;

% PLS by SIMPLS
for a = 1:A
	[~, ~, q] = svds(S'*S, 1);
	w_ast = S*q;
	t = X*w_ast;
%  	t = t - T*(T'*t);
	t = t - mean(t);
	norm_t = norm(t);
	t = t/norm_t;
	w_ast = w_ast/norm_t;
	p = X'*t;
	q = Y'*t;
	u = Y*q;
	v = p;
	v = v - V*(V'*p);
	u = u - T*(T'*u);
	v = v/norm(v);
	S = S - v*(v'*S);
% 	S = S - V*(V'*S); % double deflation to remove noise crept in
	W_ast(:, a) = w_ast;
	P(:, a) = p;
	T(:, a) = t;
	Q(:, a) = q;
	U(:, a) = u;
	V(:, a) = v;
end

switch norm_sch
	case 'standard'
		b = ones(1, A);
		U = U/diag(sqrt(diag(U'*U)));
	case 'loadings'
		norm_P = sqrt(diag(P'*P))';
		T = T*diag(norm_P);
		W_ast = W_ast*diag(norm_P);
		P = P/diag(norm_P);
		norm_Q = sqrt(diag(Q'*Q))';
		Q = Q/diag(norm_Q);
		U = U/diag(sqrt(diag(U'*U)))*diag(norm_Q);
		b = norm_Q./norm_P;
end

% % Normalisation scheme by looping
% norm_P = zeros(A, 1);
% norm_Q = zeros(A, 1);
% b = zeros(A, 1);
% for a = 1:A
%   norm_P(a) = norm(P(:, a));
%   T(:, a) = T(:, a)*norm_P(a);
%   W_ast(:, a) = W_ast(:, a)*norm_P(a);
%   P(:, a) = P(:, a)/norm_P(a);
%   norm_Q(a) = norm(Q(:, a));
%   U(:, a) = U(:, a)*(norm_Q(a)/norm(U(:, a)));
%   Q(:, a) = Q(:, a)/norm_Q(a);
%   b(a) = norm_Q(a)/norm_P(a);
% end

% Regression coefficients
B = W_ast*diag(b)*Q';

% % Loadings directed as the largest component
% [~, X_index] = max(abs(P), [], 1);
% X_colsign = sign(P(X_index + (0:V_X:(A - 1)*V_X)));
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
V_out = V;

end