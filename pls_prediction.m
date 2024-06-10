function [Y_pred_out, T_pred_out, R_new_out, E_new_out] = pls_prediction (X_new_in, Y_new_in, B_in, W_ast_in, P_in)
% PLS_PREDICTION Prediction based on PLS model
%
% Syntax:
%	[Y_pred, T_pred, R_new, E_new] = pls_prediction (X_new, Y_new, B, W_ast, P)
%
% Inputs:
%	X_new:	Array of predictors variables to apply the model to
%	Y_new:	Array of response variables to be used for error estimation (if not
%			available set it to [], corresponding errors will be set to NaNs)
%	B:		Matrix of PLS outer regression coefficients
%	W_ast:	Corrected weights matrix
%	P:		Input loadings matrix
%
% Outputs:
%	Y_pred:	Array of predicted output variables
%	T_pred:	Predicted scores matrix
%	R_new:	Prediction residuals (array of NaNs is returned in Y_new = [])
%	E_new:	Input reconstruction residuals
%
%
% NOTE
%	The data array to be used for prediction must be pre-processed in the same
%	way as the input data array used for model calibration (meaning, columns
%	must be centred on the means of the calibration dataset and scaled on the
%	variances od the calibration dataset). The array of output variables is
%	returned scaled as well.
%

%{
pls_prediction.m
Version: 1.0.0
Date: 2023-12-21
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_prediction.m
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

X_new = X_new_in;
Y_new = Y_new_in;
B = B_in;
W_ast = W_ast_in;
P = P_in;

%% PLS model prediction

% Predicted output variables
Y_pred = X_new*B;
% Predicted input scores
T_pred = X_new*W_ast;
% Reconstructued input variables
X_rec = T_pred*P';

% Check if true output variables are available
if ~isempty(Y_new)
	% Prediction residuals
	R_new = Y_new - Y_pred;
else
	% Provide NaNs
	[N, V_Y] = size(Y_pred);
	R_new = NaN(N, V_Y);
end
% Input reconstruction errors
E_new = X_new - X_rec;

%% Output assignments

Y_pred_out = Y_pred;
T_pred_out = T_pred;
R_new_out = R_new;
E_new_out = E_new;

end