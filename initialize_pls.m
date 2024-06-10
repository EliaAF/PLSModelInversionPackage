function model_out = initialize_pls ()
% INITIALISE_PLS Intilisation of the PLS model structure
%
% Syntax:
%	model = initialize_pls ()
%
% Outputs:
%	model:	Structured array of the initialized PLS model
%
%
% FIELDS OF THE STRUCTURE
%
% 	model:	main structure of the model
% 		parameters:	parameters of the PLS model
%			B:			PLS regression coefficients
% 			W_ast:		X weights (corrected)
% 			P:			X loadings
% 			b:			inner-relation regression coefficients
% 			Q:			Y loadings
% 			W:			X weights
%			sigma_sq_T:	variance of X scores for diagnostics
%			sigma_sq_U:	variance of Y scores for diagnostics
% 		prediction:	results of the application of the model
% 			T:		X scores
% 			U:		Y scores
% 			Y_pred:	predicted Y
% 			X_rec:	reconstructed X
% 			Y_rec:	reconstructed Y
% 			R:		residuals for prediction of Y
% 			E:		residuals for reconstruction of X
% 			F:		residuals for reconstruction of Y
% 		data:	data used for model calibration
% 			X:		array of X variable scaled data
% 			Y:		array of Y variable scaled data
% 			X_uns:	array of X variable unscaled data
% 			Y_uns:	array of Y variable unscaled data
% 		dimensions:	dimensions of the arrays
% 			N:		number of observations
% 			V_X:	number of X variables
% 			V_Y:	number of Y variables
% 			A:		number of latent variables
% 		scaling:	scaling parameters
% 			mu_X:		means of the X variables
% 			sigma_X:	standard deviations of the X variables
% 			mu_Y:		means of the Y variables
% 			sigma_Y:	standard deviations of the Y variables
% 		performance:	performance of the model in calibration
% 			PV_Y:			predicted variance on Y
% 			CPV_Y:			cumulative predited variance on Y
% 			EV_X:			explained variance on X
% 			CEV_X:			cumulative explained variance on X
% 			EV_Y:			explained variance on Y
% 			CEV_Y:			cumulative explained variance on Y
% 			R_sq_Y:			determination coefficient on each Y variable
% 			CR_sq_Y:		cumulative determination coefficient on each Y variable
% 			SSE:			sum of squared errors on each Y variable
% 			RMSE:			root mean squared errors on each Y variable
% 			bias:			residual biases on each Y variable
% 			SE:				standard errors on each Y variable
% 			lambda_T:		eigenvalues extracted by X scores
% 			lambda_U:		eigenvalues extracted by Y scores
% 		selection:	criteria for model selection using in-sample information
%			dof:			degrees of freedom of the model
%			R_sq_adj_Y:		adjusted determination coefficient on each Y variable
% 			IC:				Table of information criteria
% 			pIC:			Table of probabilties from information criteria
% 		diagnostics:	diagnostics on the model in calibration
% 			T_sq_T:				Hotelling's T^2 of X scores
% 			T_sq_U:				Hotelling's T^2 of Y scores
% 			SPE_Y:				squared prediction error on Y
% 			SRE_X:				squared reconstruction error on X
% 			SRE_Y:				squared reconstruction error on Y
% 			T_sq_T_con:			contribution to Hotelling's T^2 of X scores from X loadings
% 			T_sq_U_con:			contribution to Hotelling's T^2 of Y scores from Y loadings
% 			T_sq_T_con_dir:		contribution to Hotelling's T^2 of X scores from X weights
% 			SPE_Y_con:			contributions to SPE_Y
% 			SRE_X_con:			contributions to SRE_X
% 			SRE_Y_con:			contributions to SRE_Y
%			VIP:				vairable importance in projection
%			obs_lev:			leverages of observations on the model
%			R_stud:				studentized residuals for prediction of Y
%			Cook_D:				Cook's distance of observations for prediction of Y
%			lim_obs_lev:		confidence limit of leverages of observations on the model
%			lim_R_stud:			confidence limit of studentized residuals for prediction of Y
%			lim_Cook_D:			confidence limit of Cook's distance of observations for prediction of Y
%		estimates:	estimated confidence limits on diagnostics
% 			lim:				singificance level for confidence limits
%			dof:				degrees of freedom of the model
% 			lim_T_sq_T:			confidence limit of Hotelling's T^2 of X scores
% 			lim_T_sq_U:			confidence limit of Hotelling's T^2 of Y scores
% 			lim_SPE_Y:			confidence limit of squared prediction error on Y
% 			lim_SRE_X:			confidence limit of squared reconstruction error on X
% 			lim_SRE_Y:			confidence limit of squared reconstruction error on Y
% 			lim_T_sq_T_con:		confidence limit of contribution to Hotelling's T^2 of X scores from X loadings
% 			lim_T_sq_U_con:		confidence limit of contribution to Hotelling's T^2 of X scores from Y loadings
% 			lim_T_sq_T_con_dir:	confidence limit of contribution to Hotelling's T^2 of X scores from X weights
% 			lim_SPE_Y_con:		confidence limit of contribution to SPE_Y
% 			lim_SRE_X_con:		confidence limit of contribution to SRE_X
% 			lim_SRE_Y_con:		confidence limit of contribution to SRE_Y
% 			l_T:				confidence ellipsoid axes for X scores
% 			l_U:				confidence ellipsoid axes for Y scores
% 		info:	information about the model
% 			preprocessing:			kind of preprocessing used
%			algorithm:				algorithm used to calibrate the model
%			normalisation:			normalisation scheme adopted in model calibration
% 			error_based_on:			scaled on unsceld predictions and errors
%			contribution_method:	method used for computing contributions to diagnostics
%			dof_method:				method used for computing degrees of freedom
%			Tsq_lim_method:			method used for computing confidence limits on T_sq_T and T_sq_U
%			SqE_lim_method:			method used for computing confidence limits on SPE_Y, SRE_X and SRE_Y
%			con_lim_method:			method used for computing confidence limits on contributions to diagnostics
%			l_kind:					number of LVs considered for the confidence ellipse
% 			obs_names:				observation names
% 			X_var_names:			X variable names
% 			Y_var_names:			Y variable names

%{
initialize_pls.m
Version: 1.0.0
Date: 2021-09-24
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

initialize_pls.m
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

%% Structure initialization

model = struct;
	% Model paramters
	model.parameters = struct;
		model.parameters.B = [];
		model.parameters.W_ast = [];
		model.parameters.P = [];
		model.parameters.b = [];
		model.parameters.Q = [];
		model.parameters.W = [];
		model.parameters.sigma_sq_T = [];
		model.parameters.sigma_sq_U = [];
	% Scores, predictions and errors
	model.prediction = struct;
		model.prediction.T = [];
		model.prediction.U = [];
		model.prediction.Y_pred = [];
		model.prediction.X_rec = [];
		model.prediction.Y_rec = [];
		model.prediction.R = [];
		model.prediction.E = [];
		model.prediction.F = [];
	% Data used for calibration
	model.data = struct;
		model.data.X = [];
		model.data.Y = [];
		model.data.X_uns = [];
		model.data.Y_uns = [];
	% Dimensions of the model entities
	model.dimensions = struct;
		model.dimensions.N = [];
		model.dimensions.V_X = [];
		model.dimensions.V_Y = [];
		model.dimensions.A = [];
	% Sclaing applied to data
	model.scaling = struct;
		model.scaling.mu_X = [];
		model.scaling.sigma_X = [];
		model.scaling.mu_Y = [];
		model.scaling.sigma_Y = [];
	% Perfomance of the model
	model.performance = struct;
		model.performance.PV_Y = [];
		model.performance.CPV_Y = [];
		model.performance.EV_X = [];
		model.performance.CEV_X = [];
		model.performance.EV_Y = [];
		model.performance.CEV_Y = [];
		model.performance.R_sq_Y = [];
		model.performance.CR_sq_Y = [];
		model.performance.SSE = [];
		model.performance.RMSE = [];
		model.performance.bias = [];
		model.performance.SE = [];
		model.performance.lambda_T = [];
		model.performance.lambda_U = [];
	% Criteria for model selection
	model.selection = struct;
		model.selection.dof = [];
		model.selection.R_sq_adj_Y = [];
		model.selection.IC = [];
		model.selection.pIC = [];
	% Diagnostics on the model
	model.diagnostics = struct;
		model.diagnostics.T_sq_T = [];
		model.diagnostics.T_sq_U = [];
		model.diagnostics.SPE_Y = [];
		model.diagnostics.SRE_X = [];
		model.diagnostics.SRE_Y = [];
		model.diagnostics.T_sq_T_con = [];
		model.diagnostics.T_sq_U_con = [];
		model.diagnostics.T_sq_T_con_dir = [];
		model.diagnostics.SPE_Y_con = [];
		model.diagnostics.SRE_X_con = [];
		model.diagnostics.SRE_Y_con = [];
		model.diagnostics.VIP = [];
		model.diagnostics.obs_lev = [];
		model.diagnostics.R_stud = [];
		model.diagnostics.Cook_D = [];
		model.diagnostics.lim_obs_lev = [];
		model.diagnostics.lim_R_stud = [];
		model.diagnostics.lim_Cook_D = [];
	% Confidence limits
	model.estimates = struct;
		model.estimates.lim = [];
		model.estimates.dof = [];
		model.estimates.lim_T_sq_T = [];
		model.estimates.lim_T_sq_U = [];
		model.estimates.lim_SPE_Y = [];
		model.estimates.lim_SRE_X = [];
		model.estimates.lim_SRE_Y = [];
		model.estimates.lim_T_sq_T_con = [];
		model.estimates.lim_T_sq_U_con = [];
		model.estimates.lim_T_sq_T_con_dir = [];
		model.estimates.lim_SPE_Y_con = [];
		model.estimates.lim_SRE_X_con = [];
		model.estimates.lim_SRE_Y_con = [];
		model.estimates.l_T = [];
		model.estimates.l_U = [];
	% Infos on the model
	model.info = struct;
		model.info.preprocessing = '';
		model.info.algorithm = '';
		model.info.normalisation = '';
		model.info.error_based_on = '';
		model.info.diagnostics_based_on = '';
		model.info.contribution_method = '';
		model.info.dof_method = '';
		model.info.Tsq_lim_method = '';
		model.info.SqE_lim_method = '';
		model.info.con_lim_method = '';
		model.info.l_kind = '';
		model.info.obs_names = {};
		model.info.X_var_names = {};
		model.info.Y_var_names = {};

%% Output assignments

model_out = model;

end