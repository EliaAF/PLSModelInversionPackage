function t_des_cl_out = pls_inversion_null_space_uncertainty_lin (t_des_in, Q_tilde_in, alpha_in, sigma_sq_in, SSE_in, dof_in, N_in, varargin)
% PLS_INVERSION_NULL_SPACE_UNCERTAINTY_LIN Uncertainty on the null space in PLS
%	model inversion with the method by Facco et al. (2019)
%
% Syntax:
%	t_des_cl = pls_inversion_null_space_uncertainty_lin (t_des, Q_tilde, alpha, sigma_sq, SSE, dof, N)
%	t_des_cl = pls_inversion_null_space_uncertainty_lin (t_des, Q_tilde, alpha, sigma_sq, SSE, dof, N, R_Y)
%
% Inputs:
%	t_des:		Score vector estimated by PLS model inversion
%	Q_tilde:	Modified output loadings matrix of the PLS model to be inverted
%				(see notes)
%	alpha:		Significance level of the confidence limit
%	sigma_sq:	Variances of the latent variables in the calibration dataset
%	SSE:		Sum of squared errors of the output variables in calibration
%	dof:		Degrees of freedom of the PLS model
%	N:			Number of observations in the calibration dataset of the PLS
%				model
%	R_Y:		Rank on the Y_des matrix (number of independent output variables)
%
% Outputs:
%	t_des_cl:	Confidence limit of the null space at t_des
%
%
% NOTE
%	The confidence limit is provided as a vector reporting the limit for all the
%	scores from PLS model inversion. The method by Facco et al. (2015) is used.
%	If t_des has been estimated from y_des, then the uncertinaty is defined as:
%		t_des +- t_des_cl
%	Such confidence limit is then simply ``propagated'' in a linear way along the
%	direction of the null space. Note that direct inversion is used in
%	uncertainty estimation of R_Y is not provided (therefore the function assumes
%	R_Y = V_Y), while regularized direct inversion is used if R_Y is provided.

%{
pls_inversion_null_space_uncertainty_lin.m
Version: 1.0.0
Date: 2024-01-16
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_inversion_null_space_uncertainty_lin.m
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

t_des = t_des_in;
Q_tilde = Q_tilde_in;
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

% Estimate PLS ``prediction'' uncertainty at t_des
y_des_cl = pls_prediction_uncertainty(t_des, alpha, sigma_sq, SSE, dof, N);

% Determine if regularizaion is needed
if regularize
	% Regularized direct inversion
	t_des_cl = y_des_cl*(Q_tilde*reg_inv(Q_tilde'*Q_tilde, R_Y));
else
	% Direct inversion
	t_des_cl = y_des_cl*((Q_tilde*Q_tilde')\Q_tilde);
end

%% Output assignments

t_des_cl_out = t_des_cl;