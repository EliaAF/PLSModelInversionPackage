function LCS_out = pls_IBO_initialize_constraints ()
% PLS_IBO_INITIALIZE_CONSTRAINT Intilization constraint structure for PLS model
%	inversion by optimization
%
% Syntax:
%	LCS = pls_IBO_initialize_constraints ()
%
% Outputs:
%	LCS:	Structured array of the constraints for optimization
%
%
% NOTE
%	The constraint structure must be initialized using this function, then
%	populated by the user according the the constraints they wish to specify on
%	the PLS model inversion by optimization. Constraints are passed to the
%	optimization by giving LCS as input to the pls_inversion_IBO_scores function.
%	The fields of LCS are used as inputs to the pls_IBO_scores_linear_constraints
%	function. See the help of these two functions for details on how to set the
%	fields of the structure.
%
% FIELDS OF THE STRUCTURE
%
% 	LCS:	main structure of the model
% 		outputs:	constraints on output variables
%			lower_bounds:		lower bounds of outputs
% 			upper_bounds:		upper bounds of outputs
% 			equality:			equality constraints on outputs
% 			Cc_lower_bounds:	coefficients for lower bounds on linear combinations of outputs
% 			Cv_lower_bounds:	constraint values for lower bounds on linear combinations of outputs
% 			Cc_upper_bounds:	coefficients for upper bounds on linear combinations of outputs
% 			Cv_upper_bounds:	constraint values for upper bounds on linear combinations of outputs
% 			Cc_equality:		coefficients for equality constraints on linear combinations of outputs
% 			Cv_equality:		constraint values for equality constraints on linear combinations of outputs
%		inputs:		constraints on input variables
%			lower_bounds:		lower bounds of outputs
% 			upper_bounds:		upper bounds of outputs
% 			equality:			equality constraints on outputs
% 			Cc_lower_bounds:	coefficients for lower bounds on linear combinations of outputs
% 			Cv_lower_bounds:	constraint values for lower bounds on linear combinations of outputs
% 			Cc_upper_bounds:	coefficients for upper bounds on linear combinations of outputs
% 			Cv_upper_bounds:	constraint values for upper bounds on linear combinations of outputs
% 			Cc_equality:		coefficients for equality constraints on linear combinations of outputs
% 			Cv_equality:		constraint values for equality constraints on linear combinations of outputs
% 		scores:		constraints on scores (optimization variables)
%			lower_bounds:		lower bounds of scores
% 			upper_bounds:		upper bounds of scores
% 			equality:			equality constraints on scores
% 			Cc_lower_bounds:	coefficients for lower bounds on linear combinations of scores
% 			Cv_lower_bounds:	constraint values for lower bounds on linear combinations of scores
% 			Cc_upper_bounds:	coefficients for upper bounds on linear combinations of scores
% 			Cv_upper_bounds:	constraint values for upper bounds on linear combinations of scores
% 			Cc_equality:		coefficients for equality constraints on linear combinations of scores
% 			Cv_equality:		constraint values for equality constraints on linear combinations of scores

%{
pls_IBO_initialize_constraints.m
Version: 1.0.0
Date: 2023-01-12
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_IBO_initialize_constraints.m
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

LCS = struct;
	% Constraints on outputs
	LCS.outputs = struct;
		LCS.outputs.lower_bounds = [];
		LCS.outputs.upper_bounds = [];
		LCS.outputs.equality = [];
		LCS.outputs.Cc_lower_bounds = [];
		LCS.outputs.Cv_lower_bounds = [];
		LCS.outputs.Cc_upper_bounds = [];
		LCS.outputs.Cv_upper_bounds = [];
		LCS.outputs.Cc_equality = [];
		LCS.outputs.Cv_equality = [];
	% Constraints on inputs
	LCS.inputs = struct;
		LCS.inputs.lower_bounds = [];
		LCS.inputs.upper_bounds = [];
		LCS.inputs.equality = [];
		LCS.inputs.Cc_lower_bounds = [];
		LCS.inputs.Cv_lower_bounds = [];
		LCS.inputs.Cc_upper_bounds = [];
		LCS.inputs.Cv_upper_bounds = [];
		LCS.inputs.Cc_equality = [];
		LCS.inputs.Cv_equality = [];
	% Constraints on scores
	LCS.scores = struct;
		LCS.scores.lower_bounds = [];
		LCS.scores.upper_bounds = [];
		LCS.scores.equality = [];
		LCS.scores.Cc_lower_bounds = [];
		LCS.scores.Cv_lower_bounds = [];
		LCS.scores.Cc_upper_bounds = [];
		LCS.scores.Cv_upper_bounds = [];
		LCS.scores.Cc_equality = [];
		LCS.scores.Cv_equality = [];

%% Output assignments

LCS_out = LCS;