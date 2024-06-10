function [X_des_out, T_des_out, G_out] = pls_inversion_RDI (Y_des_in, Q_tilde_in, P_in, R_Y_in)
% PLS_INVERSION_RDI PLS model inversion by regularized direct inversion (RDI)
%
% Syntax:
%	X_des = pls_inversion_RDI (Y_des, Q_tilde, P, R_Y)
%	[X_des, T_des] = pls_inversion_RDI (_)
%	[X_des, T_des, G] = pls_inversion_RDI (_)
%
% Inputs:
%	Y_des:		Array of response variables to be used as targets in inversion
%	Q_tilde:	Modified output loadings matrix of the PLS model to be inverted
%				(see notes)
%	P:			Input loadings matrix of the PLS model to be inverted
%	R_Y:		Rank on the Y_des matrix (number of independent output variables)
%
% Outputs:
%	X_des:		Array of input variables estimated by model inversion
%	T_des:		Scores matrix estimated by model inversion
%	G:			Basis matrix of the null space of the PLS model
%
%
% NOTE
%	The data array to be used as output target must be pre-processed in the same
%	way as the output data array used for model calibration (meaning, columns
%	must be centred on the means of the calibration dataset and scaled on the
%	variances of the calibration dataset). The array of input variables is
%	returned scaled as well.
%
% NOTE
%	The inversion method implemented in this function is based on the simplified
%	PLS model commonly used in the literature:
%		T = X*W_ast
%		X_est = T*P'
%		Y_est = T*Q_tilde'
%	where matrix Q_tilde is obtained by the output loading matrix (Q) and inner
%	regression coefficients (b) of the general PLS model as:
%		Q_tilde = Q*diag(b)
%	The function expects matrix Q_tilde to be provided as input.
%
% NOTE
%	The RDI method assumes that output variables can be correlated, thus only R_Y
%	variables out of V_Y variables are independent.If there are more latent
%	variables than independent output variables, a null space exists and the
%	inversion has infinite solutions. The function returns only one "optimal"
%	solution (in the sense of the Frobenious Norm) and matrix G, the columns of
%	which form a basis of the null space. If the number of latent variables is
%	less of equal to the number of independent output variables, the solution to
%	model inversion in unique and matrix G is returned as [].

%{
pls_inversion_RDI.m
Version: 1.0.0
Date: 2023-12-21
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

pls_inversion_RDI.m
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

Y_des = Y_des_in;
Q_tilde = Q_tilde_in;
P = P_in;
R_Y = R_Y_in;

%% PLS model inversion

% Number of latent variables
A = size(Q_tilde, 2);

% Invert model
T_des = Y_des*(Q_tilde*reg_inv(Q_tilde'*Q_tilde, R_Y));

% Determine case for null space
if A <= R_Y
	% No null space
	G = [];
else
	% Yes null space
	[G, ~, ~] = svd(Q_tilde');
	G(:, 1:R_Y) = [];
end
% Compute input variables
X_des = T_des*P';

%% Output assignments

X_des_out = X_des;
if nargout > 1
	T_des_out = T_des;
end
if nargout > 2
	G_out = G;
end

end