function [X, Y] = ellipse (x0, y0, a, b)
% ELLIPSE calculates vectors to draw a 2D ellipse
%
% Syntax:
%	[X, Y] = ellipse (x0, y0, a, b)
%
% Inputs:
%	x0:		Abscissa coordinate of the center
%	y0:		Ordinate coordinate of the center
%	b:		Major semiaxis
%	a:		Minor semiaxis
%
% Outputs:
%	X:		Abscissas vector for plotting
%	Y:		Ordinates vector for plotting

%{
ellipse.m
Version: 1.0.0
Date: 2021-04-02
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

ellipse.m
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

% Parameter vector in degrees
alpha = linspace(0, 360, 360)';
% Parameter vector in radians
alpha = alpha.*(pi/180);
% Sine and cosine of the parameter
sinalpha = sin(alpha);
cosalpha = cos(alpha);
% Vectors to be plotted
X = x0 + a*cosalpha;
Y = y0 + b*sinalpha;
end