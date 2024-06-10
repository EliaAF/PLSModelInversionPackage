function T_null_cl_out = pls_inversion_null_space_uncertainty_nonlin (T_null_in, Q_tilde_in, y_des_in, alpha_in, sigma_sq_in, SSE_in, dof_in, N_in, varargin)
% PLS_INVERSION_NULL_SPACE_UNCERTAINTY_NONLIN Uncertainty on the null space in
%	PLS model inversion with the method by Palací-López et al. (2019)
%
% Syntax:
%	T_null_cl = pls_inversion_null_space_uncertainty_nonlin (T_null, Q_tilde, y_des, alpha, sigma_sq, SSE, dof, N)
%	T_null_cl = pls_inversion_null_space_uncertainty_nonlin (T_null, Q_tilde, y_des, alpha, sigma_sq, SSE, dof, N, R_Y)
%
% Inputs:
%	T_null:		Matrix of scores alonf the null space in PLS model inversion
%	Q_tilde:	Modified output loadings matrix of the PLS model to be inverted
%	y_des:		Target used in PLS model inversion
%	alpha:		Significance level of the confidence limit
%	sigma_sq:	Variances of the latent variables in the calibration dataset
%	SSE:		Sum of squared errors of the output variables in calibration
%	dof:		Degrees of freedom of the PLS model
%	N:			Number of observations in the calibration dataset of the PLS
%				model
%	R_Y:		Rank on the Y_des matrix (number of independent output variables)
%
% Outputs:
%	T_null_cl:	Array of confidence limit of the null space at each point in
%				T_null
%
%
% NOTE
%	The output target vector must be pre-processed in the same way as the output
%	data array used for model calibration (meaning, columns must be centred on
%	the means of the calibration dataset and scaled on the variances of the
%	calibration dataset).
%
% NOTE
%	This function requires an array of points along the null space as input.
%	Every row in T_null is regarded as an ``observation'' along the null space.
%	The confidence limit is provided as an array reporting the limit for all the
%	points in T_null. The method by Palací-López et al. (2019) is used. The
%	uncertinaty is defined as:
%		T_null +- T_null_cl
%	Note that direct inversion is used in uncertainty estimation of R_Y is not
%	provided (therefore the function assumes R_Y = V_Y), while regularized direct
%	inversion is used if R_Y is provided.

%{
pls_inversion_null_space_uncertainty_nonlin.m
Version: 1.0.0
Date: 2024-01-16
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_inversion_null_space_uncertainty_nonlin.m
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

T_null = T_null_in;
Q_tilde = Q_tilde_in;
y_des = y_des_in;
alpha = alpha_in;
sigma_sq = sigma_sq_in;
SSE = SSE_in;
dof = dof_in;
N = N_in;

%% Optional arguments development

% Optionals initialization
regularize = false;

% Development cycle
if ~isempty(varargin)
	regularize = true;
	R_Y = varargin{1};
end

%% Null space uncertainty

% Residual model error
R_null = y_des - T_null*Q_tilde';
% Determine if regularizaion is needed
if regularize
	% Regularized direct inversion
	T_tilde = T_null + R_null*(Q_tilde*reg_inv(Q_tilde'*Q_tilde, R_Y));
else
	% Direct inversion
	T_tilde = T_null + R_null*((Q_tilde*Q_tilde')\Q_tilde);
end
% Estimate PLS ``prediction'' uncertainty at T_tilde
Y_des_cl = pls_prediction_uncertainty(T_tilde, alpha, sigma_sq, SSE, dof, N);

% Determine if regularizaion is needed
if regularize
	% Regularized direct inversion
	T_null_cl = Y_des_cl*(Q_tilde*reg_inv(Q_tilde'*Q_tilde, R_Y));
else
	% Direct inversion
	T_null_cl = Y_des_cl*((Q_tilde*Q_tilde')\Q_tilde);
end

%% Output assignments

T_null_cl_out = T_null_cl;