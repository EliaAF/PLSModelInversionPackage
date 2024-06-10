function [IC_out, pIC_out] = information_criteria (SSE_in, dof_in, N_in)
% INFORMATION_CRITERIA Model selection based on information criteria
%
% Syntax:
%	IC = information_criteria (SSE, dof, N)
%	[IC, pIC] = information_criteria (SSE, dof, N)
%
% Inputs:
%	SSE:	Column vector of sum of squared errors of models to be compared
%	dof:	Column vector of degrees of freedom of models to be compared
%	N:		Number of observations in the original dataset used for model fitting
%
% Outputs:
%	IC:		Table of values of information criteria
%	pIC:	Table of probabilities of selecting the best model according to
%			infomation criteria
%
%
% NOTE
%	Information criteria considered in the IC and pIC tables
%		FPE:		Akaike's Final Prediction Error
%		Cp:			Mallows' Cp
%		GMC:		Geweke-Meese Criterion
%		AIC:		Akaike Information Criterion
%		BIC:		Bayesian Information Criterion
%		HCQ:		Hannah-Quinn Criterion
%		AICc:		Akaike Information Criterion corrected (small sample)
%		BICc:		Bayesian Information Criterion corrected (small sample)
%		HCQc:		Hannah-Quinn Criterion corrected (small sample)
%	Noote that also SSE and dof are collected into the IC table
%
% NOTE
%	Probabilities are computed on offsets of IC values on the minimum among the
%	ones to be compared, then assumed to follow an exponential distribution

%{
information_criteria.m
Version: 1.0.0
Date: 2022-01-28
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

information_criteria.m
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

SSE = SSE_in;
dof = dof_in;
N = N_in;

%% Utilities

% Information criteria to be considered
IC_names = {
	'SSE'
	'dof'
	'FPE'
	'Cp'
	'GMC'
	'AIC'
	'BIC'
	'HQC'
	'AICc'
	'BICc'
	'HQCc'
};
p_IC_names = cellfun(@(b) ['p' b], IC_names(3:end), 'UniformOutput', false);

% Initialise table of information criteria
IC = table('Size', [size(SSE, 1), size(IC_names, 1)], 'VariableTypes', repmat({'double'}, size(IC_names, 1), 1));
IC.Properties.VariableNames = IC_names;

% % Initialise table of probabilties
% pIC = table('Size', [size(SSE, 1), size(p_IC_names, 1)], 'VariableTypes', repmat({'double'}, size(p_IC_names, 1), 1));
% pIC.Properties.VariableNames = p_IC_names;

% Save inputs into the IC table
IC.SSE = SSE;
IC.dof = dof;

%% Information criteria

% MSE in calibration (prediction error on Y)
MSE = SSE/N;
% Unbiased MSE in calibration (prediction error on Y)
MSE_unb = SSE./(N - dof);

% % Adjusted determination coefficient on each response variable
% IC.R_sq_adj = 1 - ((N - 1)./(N - dof)).*(1 - R_sq_Y);

% Akaike's Final Prediction Error
IC.FPE = MSE.*(N + dof)./(N - dof);
% Mallows' Cp
IC.Cp = SSE./MSE_unb + 2*dof;
% Geweke-Meese Criterion
IC.GMC = SSE./MSE_unb + dof*log(N);

% Akaike Information Criterion
IC.AIC = N*log(MSE) + 2*dof;
% Bayesian Information Criterion
IC.BIC = N*log(MSE) + dof*log(N);
% Hannah-Quinn Criterion
IC.HQC = N*log(MSE) + 2*dof*log(log(N));

% Corrector function
cor_fun = (dof + 1)./(N - (dof + 1));

% Akaike Information Criterion corrected
IC.AICc = IC.AIC + 2*dof.*cor_fun;
% Bayesian Information Criterion corrected
IC.BICc = IC.BIC + dof*log(N).*cor_fun;
% Hannan-Quinn Criterion corrected
IC.HQCc = IC.HQC + 2*dof*log(log(N)).*cor_fun;

% Probability distributions from model selection criteria
pIC = varfun(@prob_dist_IC, IC(:, 3:end));
pIC.Properties.VariableNames = p_IC_names;

%% Output assignments

IC_out = IC;
pIC_out = pIC;
