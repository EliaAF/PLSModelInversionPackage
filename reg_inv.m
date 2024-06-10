function A_inv = reg_inv (A, R)
% REG_INV Regularized inverse of a matrix computed by the Singular Value
%	Decomposition (SVD) algorithm with regularaziation
%
% Syntax:
%	A_inv = reg_inv (A, R)
%
% Inputs:
%	A:		Matrix to be inverted
%	R:		Number of singular values to be used for regularised inversion
%			(should match the true rank of A)
%
% Outputs:
%	A_inv:	Generalised inverse matrix (with regularisation)
%
%
% NOTE
%	The generalized inverse of a matrix, also called pseudo-inverse or
%	Moore-Penrose inverse, is an extension of inverse matrices to non-square
%	matrices. This function computes such an inverse using the SVD algorithm,
%	which is also exploited for regularizing the inversion. Namely, only the
%	R larger singlular values are used in the inversion (if R is the rank of A,
%	such singular values are the non-null ones). Is the specified R is larger
%	than the minimum dimension of A, R is set to A automatically.

%{
reg_inv.m
Version: 1.0.0
Date: 2021-09-21
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

reg_inv.m
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

% Check if R is feasible
R_max = min(size(A));
if R > R_max
	R = R_max;
end

% Singular value decomposition of A
[U, S, V] = svd(A);

% Regularization
V = V(:, 1:R);
S = S(1:R, 1:R);
U = U(:, 1:R);

% Regularised inverse computation
A_inv = (V/S)*U';
end