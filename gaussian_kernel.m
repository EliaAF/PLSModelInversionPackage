function p = gaussian_kernel (x, x_n, h)
% GAUSSIAN_KERNEL Univariate Gaussian kernel for density estimation
%
% Syntax:
%	p = gaussian_kernel (x, x_n, h)
%
% Inputs:
%	x:		Evaluation point
%	x_n:	Central point of the kernel
%	h:		Kernel bandwidth
%
% Outputs:
%	p:		Value of the probability density function of the kernel at x
%
%
% NOTE
%	The function can be used to compute the density function of several
%	evaluation points for kernels at severl central points with same bandwidth,
%	which is particularly helpful for kernel density estimation. Provide the
%	evaluation point at a column vector and the central points as a row vector:
%	p is returned as a matrix containing the PDF of kenels centred on all the
%	points provided (rows are evaluation points, columns are central points).

%{
gaussian_kernel.m
Version: 1.0.0
Date: 2023-12-20
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

gaussian_kernel.m
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

% Compute the probability density function of the Gaussin Kernel
p = exp(-0.5*((x - x_n).^2)/h^2)/sqrt(2*pi);
end