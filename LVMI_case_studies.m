% Primary script for PLS model calibration and inversion

%{
LVMI_case_studies.m
Version: 1.0.0
Date: 2024-02-20
Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

This file is part of the LVM codebase developed by Elia Arnese Feffin, and is
covered by the GNU General Public License version 3

# GNU General Public License version 3 (GPL-3.0) --------------------------------

LVMI_case_studies.m
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

%% Instructions and information

%{
This script guides the user through the workflow for PLS model calibration,
selection, and inversion. Three methods for model inversion (direct inversion,
regularized direct inversion, and inversion by optimization) are carried out and
compared. Inputs to be reviewed and/or set by the user are distributed throughout
this script. Descriptions of the actions to be performed to use
this script are marked in the code as:
	%---> DESCRIPTION OF ACTION <---%
in the immediate proximity of the lines of code to be changed. Additional
information are provided under the instruction line, if appropriate. As inputs
are not gathered in a single section and are spread in the code, it is warmly
recommended to execute the code section by section, reviewing all the
instructions in doing so.

The general structure of the script is as follows.
	- Section ``Data section'': dataset loading and output target assignment.
	- Section ``Analysis parameters'': main settings of the PLS model selection
	  and calibration procedures.
	- Section ``Model assessment and selection'': PLS model selection routines.
	- Section ``Model calibration'': PLS model calibration routines.
	- Section ``Inversion: common settings'': main settings of the PLS model
	  inversion procedures.
	- Section ``Direct inversion (DI)'': PLS model inversion routines for direct
	  inversion.
	- Section ``Redularized direct inversion (RDI)'': PLS model inversion
	  routines for regularized direct inversion.
	- Section ``Inversion by optimization (IbO)'': PLS model inversion routines
		for inversion by optimization, including settings of the optimization
		algorithms (e.g., weights for the soft constraints, hard constraints to
		be stated, etc.).
	- Section ``Null space'': determination of the null space in PLS model
	  inversion, if it exists.
	- Section ``Plots'': Plots for results.
Throughout this script, the columns of data matrices X and Y are assumed to be
the input and output variables, respectively, while their rows are the input and
output observations, respectively. The number of input variables is V_X, the
number of output variables is V_Y, and the number of latent variables in the PLS
model is A. Vectors (observations) of all variables are assumed to be row
vectors: considering, for example, a vector of input variables x, then
size(x) = [1, V_X]. The values of constraints in IbO are to be passed as row
vectors as well.

This script offers two pre-implemented case studies:
	- The ``example data'' provided by Palací-López et al. (2019), Chemometr.
	  Intell. Lab. Syst. 194, 103848.
	- The ``bioreactor case study'' discussed by Arnese-Feffin et al. (2022),
	  Chemometr. Intell. Lab. Syst. 227, 104610.
Additionally, this script allows the user to provide custom data for PLS model
calibration and inversion. In this case, please carefully review the instructions
spread throughout the script. Note that all the case studies are gathered in a
single script for ease of execution: case-specific settings are contained in
switch-case statements throughout the code.

Bear in mind that this script is to be interpreted as an example rather than a
general interface, especially regarding the ``Plots'' section. Consider the
following limitations of this script.
	- The pre-implemented case studies should be executed using A = 2 LVs; errors
	  will be obtained in the ``Plots'' section otherwise.
	- The code in the ``Plots'' section is meant to handle only the case of a
	  1-dimensional null space.
	- No plots are produced if custom data are used.
	- No check on the feasibility of constraints in IbO is performed.
%}


%% Data selection

%---> UN-COMMENT THE DATASET TO BE USED <---%
% If custom data are used, review carefully all the instructions and set all the
% required inputs.

% Selector for dataset
datasel = 'palaci-lopez';	% Example data by Palací-López et al. (2019), Chemometr. Intell. Lab. Syst. 194, 103848
% datasel = 'bioreactor';	% Bioreactor case-study data by Arnese-Feffin et al. (2022), Chemometr. Intell. Lab. Syst. 227, 104610
% datasel = 'custom';		% Custom dataset provided by the user (as .mat file)

% Select dataset
switch datasel
	case 'palaci-lopez'
		% Load data
		load dataset_palacilopez2019_ex.mat
		% Numbers of observations and of input variables
		[N, V_X] = size(X);
		% Number of output variables
		V_Y = size(Y, 2);

		%---> UNCOMMENT THE TARGET TO BE USED <---%
		% In this case study, the target is computed by first specifying two
		% independent input variables by the row vector idx, then transforming it
		% to obtain all five inputs (x_des), and finally computing the output
		% target (y_des) by the true linear model.

		% Assign independent variables for quality target
		idx = [36.27, 10.80];	% Reasonable target
% 		idx = [0, 1];			% Unreasonably low target to check hard contraints
% 		idx = [66.27, 22.80];	% Unreasonably high target to check hard contraints

		% Form input vector for the target
		x_des = [idx(1), idx(2), idx(1)^2, idx(2)^2, idx(1)*idx(2)];
		% Compute output target
		y_des = x_des*[4.3; 0.022; -0.0064; 1.1; -0.12] - 21;

	case 'bioreactor'
		% Load data
		load dataset_batchcellculture.mat
		% Numbers of observations and of input variables
		[N, V_X] = size(X);
		% Number of output variables
		V_Y = size(Y, 2);
		% Index of the observation to be used as output target
		idx = 13;
		% Assign output target
		y_des = Y(idx, :);
		% Assign true input variables
		x_des = X(idx, :);
		% Remove target observation from the dataset
		X(idx, :) = [];
		Y(idx, :) = [];
		% Adjust the number of observations
		N = N - 1;

	case 'custom'

		%---> PROVIDE NAME OF THE CUSTOM DATASET TO BE LOADED <---%
		% The .mat file to be used as custom dataset should contain (at least)
		% matrices for the input and ouptut variables named X and Y,
		% respectively. Columns of the matrices represent variables, while rows
		% represent observations.

		% Load data
		load ''

		% Numbers of observations and of input variables
		[N, V_X] = size(X);
		% Number of output variables
		V_Y = size(Y, 2);

		%---> PROVIDE OUTPUT TARGET AS ROW VECTOR <---%
		% The target output is to be provided as a variable named y_des, a row
		% vector such that size(y_des) = [1, V_Y] whose columns constitute the
		% targets on single output variables. If some output variables do not
		% have a target value, set the corresponding component of y_des to NaN
		% and adjust the vector y_free in section
		% ``Inversion by optimization (IbO)''.

		% Assign output target
		y_des = zeros(1, V_Y);
		
		%---> PROVIDE TRUE INPUTS AS ROW VECTOR <---%
		% The true input vectore corresponding to y_des can be provided as a
		% variable named x_des, a row vector such that size(x_des) = [1, V_X].
		% Set ``x_des = NaN(1, V_X)'' if the true input vector is unknown.

		% Assign true input variables
		x_des = zeros(1, V_X);
end

%% Analysis parameters

%---> UNCOMMENT THE MODEL SELECTION METHOD TO BE USED <---%
% For the scope of the present script, it is recommended to manually assign the
% number of latent variables. Note that the plots implemented in the ``Plots''
% section of this script can handle cases 2 LVs for both
% datasel == 'palaci-lopez' and datasel == 'bioreactor'.

% Selector for the number of LVs
A_selector = {
% 	'cross-validation'			% Cross-validation
% 	'cross-validation_oster'	% Cross-validation with one-stadard-error rule
	'manual'					% Manual assignment
};
A_selector = char(A_selector);

%---> REVIEW OR FILL-IN THE FOLLOWING SETTINGS FOR CROSS-VALIDATION <---%
% See the help of the function ``cross_validate_pls.m'' for details

% Number of latent-variables to cross-validate
A_CV = min(V_X, N) - 1;
% Cross-validation splitting method
CV_method = 'random_subsets';
% Number of cross-validation splits
G_obs = 5;
% Number of reptitions of random_subsets corss-validation
N_rep = 10;
% Parsimony factor for one-standard-error rule
pfac = 1;

%---> ASSIGN THE DESIRED NUMBER OF LVs FOR MANUAL ASSIGNMENT <---%

% Manual assignment of number of latent variables
A_man = 2;

%---> ASSIGN THE CONFIDENCE LIMIT FOR DIAGNOSTICS <---%

% Confidence limit for statistics
lim = 0.95;

%% Model assessment and selection

%---> REVIEW, CHANGE, OR EXTEND SETTINGS FOR PLS CROSS-VALIDATION <---%
% See the help of the function ``cross_validate_pls.m'' for details

% Model selection based on cross-validation
model_CV = cross_validate_pls(X, Y, A_CV, CV_method, G_obs, N_rep,...
	'Algorithm', 'simpls',...
	'Aggregate', false,...
	'ParFac', pfac,...
	'Normalise', 'loadings',...
	'ConLimCV', lim,...
	'DOFMethodCV', 'naive',...
	'LimMethodCV', 'norm'...
);

% Overall optimal number of LVs
[~, A_opt_CV] = min(mean(model_CV.statistics.MSE_split.pop, [2, 3, 4]), [], 1);
% Optimal number of LVs for each Y variables
[~, A_opt_CV_V_Y] = min(mean(model_CV.statistics.MSE_split.pop, [3, 4]), [], 1);

% Overall optimal number of LVs with one-standard-error-rule
mean_MSE = mean(model_CV.statistics.MSE_split.pop, [2, 3, 4]);
std_MSE = std(model_CV.statistics.MSE_split.pop, 0, [2, 3, 4]);
min_MSE = min(mean_MSE, [], 1);
idx = mean_MSE <= min_MSE + (pfac*std_MSE)/sqrt(G_obs*N_rep);
A_opt_CV_oster = find(idx, 1);
% Optimal number of LVs for each Y variables with one-standard-error-rule
A_opt_CV_oster_V_Y = model_CV.results.A_sug;

%% Model calibration

% Assignment of the number of principal components according to model selection
switch A_selector
	case 'cross-validation'
		A = A_opt_CV;
	case 'cross-validation_oster'
		A = A_opt_CV_oster;
	case 'manual'
		A = A_man;
end

%---> REVIEW, CHANGE, OR EXTEND SETTINGS FOR PLS CALIBRATION <---%
% See the help of the function ``build_pls.m'' for details

% Calibrate PLS model
model = build_pls(X, Y, A,...
	'Algorithm', 'simpls',...
	'Normalise', 'loadings',...
	'ErrBasedOn', 'scaled',...
	'DOFMethod', 'naive'...
);

% Get entities
B = model.parameters.B;
W_ast = model.parameters.W_ast;
Q = model.parameters.Q;
P = model.parameters.P;
b = model.parameters.b;
sigma_sq = model.parameters.sigma_sq_T;

T = model.prediction.T;

EV_X = model.performance.EV_X;
PV_Y = model.performance.PV_Y;
SSE = model.performance.SSE(A, :);
CR_sq_Y = model.performance.CR_sq_Y(A, :);

T_sq = model.diagnostics.T_sq_T;
SRE = model.diagnostics.SRE_X;

lim_T_sq = model.estimates.lim_T_sq_T;
lim_SRE = model.estimates.lim_SRE_X;
l_T = model.estimates.l_T;
dof = model.estimates.dof;

%% Inversion: common settings

% Scale quality target and true input
y_des_s = scale_by(y_des, model.scaling.mu_Y, model.scaling.sigma_Y);
x_des_s = scale_by(x_des, model.scaling.mu_X, model.scaling.sigma_X);

% Compute true scores of target
t_des = x_des_s*W_ast;
% Matrix Q_tilde
Q_tilde = Q*diag(b);
% Significance level from confidence limit
alpha = 1 - lim;

%% Direct inversion (DI)

% PLS model direct inversion
[x_des_DI, t_des_DI] = pls_inversion_DI(y_des_s, Q_tilde, P);
% Predicted quality
y_des_DI = rescale_by(t_des_DI*Q_tilde', model.scaling.mu_Y, model.scaling.sigma_Y);
% Diagnostics of the solution
T_sq_DI = T_sq_statistic(t_des_DI, sigma_sq);
SRE_DI = Q_statistic(x_des_DI*(eye(V_X) - W_ast*P'));
% Re-scale designed input variables
x_des_DI = rescale_by(x_des_DI, model.scaling.mu_X, model.scaling.sigma_X);

%% Regularized direct inversion (RDI)

% Rank for RDI
switch datasel
	case 'palaci-lopez'
		R_Y = 1;
	case 'bioreactor'
		R_Y = 1;
	case 'custom'
		R_Y = 0;	%---> ASSIGN THE NUMBER OF INDEPENDENT OUTPUT VARIABLES <---%
end

% PLS model regularized direct inversion
[x_des_RDI, t_des_RDI] = pls_inversion_RDI(y_des_s, Q_tilde, P, R_Y);
% Predicted quality
y_des_RDI = rescale_by(t_des_RDI*Q_tilde', model.scaling.mu_Y, model.scaling.sigma_Y);
% Diagnostics of the solution
T_sq_RDI = T_sq_statistic(t_des_RDI, sigma_sq);
SRE_RDI = Q_statistic(x_des_RDI*(eye(V_X) - W_ast*P'));
% Re-scale designed input variables
x_des_RDI = rescale_by(x_des_RDI, model.scaling.mu_X, model.scaling.sigma_X);

%% Inversion by optimization (IbO)

%---> REVIEW OR MODIFY SETTINGS FOR INVERSION BY OPTIMIZATION BELOW <---%
% See the function ``pls_inversion_IBO.m'' for details

%---> SPECIFY OUTPUT VARIABLES WITHOUT A TARGET VALUE <---%
% Set the i-th component of ``y_free'' to ``true'' if the i-th output variable
% does not have a target value.

% Free outout variables
y_free = false(1, V_Y);
% y_free(1) = true;

% Assign target for IbO
y_targ_IbO = y_des_s;
% Free some output variables
y_targ_IbO(y_free) = NaN;

%---> SELECT OPTIMIZATION DOMAIN <---%
% Set ``nonzero_SRE = false'' to optimize in the space of latent variables or
%``nonzero_SRE = true'' to optimize in the space of input variable.

% Allow non-null SRE of the inversion solution
nonzero_SRE = false;

%---> REVIEW OR MODIFY WEIGHTS OF THE OBJECTIVE FUNCTION <---%
% Gamma_1 must be provided as a row vector such that size(Gamma_1) = [1, V_Y],
% where the components of Gamma_1 are the weights of individual output variables
% on the objective function of the IbO problem. gamma_2 and gamma_3 are scalar
% values acting as weights for the soft constraints on T^2 and Q (here denoted as
% squared-reconstruction error, SRE) of the IbO solution, respectively.

% Weight vector for the targets on output variables
Gamma_1 = CR_sq_Y;		% Set ``Gamma_1 = ones(1, V_Y) for the unweighted case
% Null weights of free outputs
Gamma_1 = Gamma_1.*(~y_free);
% Weight for soft constraint on T^2
gamma_2 = 1/lim_T_sq;	% Set ``gamma_2 = 0'' to remove the soft constraint on T^2
% Weight for soft constraint on SRE
gamma_3 = 1/lim_SRE;	% Set ``gamma_3 = 0'' to remove the soft constraint on SRE

%---> REVIEW OR MODIFY LINEAR CONSTRAINTS ON OUTPUS, INPUTS, AND SCORES <---%
% Constraints are to be provided as row vectors of approriate numbers of
% components. For conveninence of the user, contraints values for outputs,
% inputs, and latent variables have been initialized as NaN vectors of the
% appropriate dimensions for the two pre-implemented case studies. Simply replace
% the NaN in the component to be constrained with the constraint value to be set.
% Note that the NaN vectors to be modified are inside functions ``scale_by'' for
% the constraint values to be adapated to the PLS model: in this script
% constraints are to be stated as unscaled variables and the scaling is done
% automatically. Constraints constraints on linear combinations of the variables
% can be set as well. See the helps of the functions
% ``pls_IBO_initialize_constraints.m'' and ``pls_IBO_linear_constraints.m'' for
% additional details.

% Initialize structure for linear constraints
LCS = pls_IBO_initialize_constraints();
% State linear constraints for the optimization problem
switch datasel
	case 'palaci-lopez'
		LCS.outputs.lower_bounds = scale_by(NaN, model.scaling.mu_Y, model.scaling.sigma_Y);
		LCS.outputs.upper_bounds = scale_by(NaN, model.scaling.mu_Y, model.scaling.sigma_Y);
		LCS.outputs.equality = scale_by(NaN, model.scaling.mu_Y, model.scaling.sigma_Y);
		LCS.inputs.lower_bounds = scale_by([NaN, NaN, NaN, NaN, NaN], model.scaling.mu_X, model.scaling.sigma_X);
		LCS.inputs.upper_bounds = scale_by([NaN, NaN, NaN, NaN, NaN], model.scaling.mu_X, model.scaling.sigma_X);
		LCS.inputs.equality = scale_by([NaN, NaN, NaN, NaN, NaN], model.scaling.mu_X, model.scaling.sigma_X);
		LCS.scores.lower_bounds = [NaN, NaN];
		LCS.scores.upper_bounds = [NaN, NaN];
	case 'bioreactor'
		LCS.outputs.lower_bounds = scale_by([NaN, NaN], model.scaling.mu_Y, model.scaling.sigma_Y);
		LCS.outputs.upper_bounds = scale_by([NaN, NaN], model.scaling.mu_Y, model.scaling.sigma_Y);
		LCS.outputs.equality = scale_by([NaN, NaN], model.scaling.mu_Y, model.scaling.sigma_Y);
		LCS.inputs.lower_bounds = scale_by([NaN, NaN, NaN, NaN], model.scaling.mu_X, model.scaling.sigma_X);
		LCS.inputs.upper_bounds = scale_by([NaN, NaN, NaN, NaN], model.scaling.mu_X, model.scaling.sigma_X);
		LCS.inputs.equality = scale_by([NaN, NaN, NaN, NaN], model.scaling.mu_X, model.scaling.sigma_X);
		LCS.scores.lower_bounds = [NaN, NaN];
		LCS.scores.upper_bounds = [NaN, NaN];
	case 'custom'
		%---> ASSIGN LINEAR CONSTRAINTS FOR THE CUSTOM DATA CASE <---%
		% Constraints are to be provided as row vectors of approriate numbers of
		% components. See the helps of the functions
		% ``pls_IBO_initialize_constraints.m'' and
		% ``pls_IBO_linear_constraints.m'' for details.
end

%---> REVIEW OR MODIFY LINEAR CONSTRAINTS ON THE NULL SPACE (IF ANY) <---%
% Additional linear constraints can be stated to snap the IbO solution to the
% null space (if it exists) or to bound it within its confidence region estimated
% by the method proposed by Facco et al. (2015), Ind. Eng. Chem. Res., 54,
% 5128−5138. See the help of the function ``pls_inversion_IBO.m'' for additional
% details.

% Snap the solution to the null space (if any)
NS_eq = false;
% Bound the solution within the null space linear confidence region
NS_lin_bnd = false;

%---> REVIEW OR MODIFY NONLINEAR CONSTRAINTS <---%
% Nonlinear constraints can be state to bound the IbO solution within the PLS
% model validity region definied by the confidence limits of the T^2 and Q (SRE)
% statistics. Additionally, a nonlinear constrain can be stated to bound the IbO
% solution within the confidence region of the null space (if it exists)
% estimated by the method proposed by Palací-López et al. (2019), Chemometr.
% Intell. Lab. Syst. 194, 103848. See the help of the function
% ``pls_inversion_IBO.m'' for additional details.

% Bound the solution withing the confidence region of T^2
T_sq_bnd = false;
% Bound the solution withing the confidence region of SRE
SRE_bnd = false;
% Bound the solution within the null space nonlinear confidence region
NS_nl_bnd = false;

% PLS model inversion by optimization
[x_des_IbO, t_des_IbO] = pls_inversion_IBO(...
	y_targ_IbO,...
	Q_tilde,...
	P,...
	W_ast,...
	sigma_sq,...
	Gamma_1,...
	gamma_2,...
	gamma_3,...
	nonzero_SRE,...
	LCS,...
	NS_lin_bnd,...
	NS_lin_bnd,...
	alpha,...
	SSE,...
	dof,...
	N,...
	[],...	%---> OPTIONS TO THE QUADRATIC PROGRAMMING SOLVER <---%
	T_sq_bnd,...
	lim_T_sq,...
	SRE_bnd,...
	lim_SRE,...
	NS_nl_bnd,...
	[]...	%---> OPTIONS TO THE NONLINEAR PROGRAMMING SOLVER <---%
);
% Predicted quality
y_des_IbO = rescale_by(t_des_IbO*Q_tilde', model.scaling.mu_Y, model.scaling.sigma_Y);
% Diagnostics of the solution
T_sq_IbO = T_sq_statistic(t_des_IbO, sigma_sq);
SRE_IbO = Q_statistic(x_des_IbO*(eye(V_X) - W_ast*P'));
% Re-scale designed input variables
x_des_IbO = rescale_by(x_des_IbO, model.scaling.mu_X, model.scaling.sigma_X);

%% Null space

% NOTE: in case some output variables are left free, the concept of null space
% becomes quite triky. In this implementation of PLS model inversion, constraints
% on the null space are allowed even if not all outputs have a target value. A
% ``reduced'' output loading matrix is used in this case, where rows
% corresponding to the free outputs are neglected.

% Check if null space exists
if A > V_Y - sum(y_free)
	% Reduced target vector
	y_targ_red_IbO = y_targ_IbO;
	y_targ_red_IbO(y_free) = 0;
	% Reduced loading matrix
	Q_tilde_red = Q_tilde;
	Q_tilde_red(y_free, :) = 0;
	% Compute reduced DI solution and null space
	t_red_DI = y_targ_red_IbO*(Q_tilde_red*reg_inv(Q_tilde_red'*Q_tilde_red, V_Y - sum(y_free)));
	[G_red, ~, ~] = svd(Q_tilde_red');
	G_red(:, 1:(V_Y - sum(y_free))) = [];

	% Null space uncertainty, method by Facco et al. (2015), Ind. Eng. Chem. Res. 54, 5128-5138
	t_red_DI_cl = pls_inversion_null_space_uncertainty_lin(t_red_DI, Q_tilde_red, alpha, sigma_sq, SSE, dof, N, V_Y - sum(y_free));
	% Confidence limits of the null space
	t_red_DI_lcl = t_red_DI - t_red_DI_cl;
	t_red_DI_ucl = t_red_DI + t_red_DI_cl;
	
	% NOTE: The following code estimates the null space based on the method by
	% Palací-López et al. (2019), Chemometr. Intell. Lab. Syst. 194, 103848,
	% which is numerically rendered by a discretization of the null space.
	% However, the discretization method implemented in this script is
	% appropriate only when the null space 1-dimensional, therefore up to 2 LVs
	% if datasel == 'palaci-lopez' and 3 LVs if datasel == 'bioreactor'. While
	% this is enough for simple visualization of the constraints, as for the
	% scope of this example script, it is by no mean a general method.

	% Extremes for the discretization of the null space
	switch datasel
		case 'palaci-lopez'
			extr = [-5, 5];
		case 'bioreactor'
			extr = [-3, 3];
		case 'custom'
			%---> STATE EXTREMES FOR THE DISCRETIZATION OF THE NULL SPACE <---%
			% The 1-dimensional null space is discretized in N_NS point (see few
			% lines below) with extremes defined by the vector extr: the first
			% component is the lower bound, while the second compone t is the
			% upper bound. Note that the comoents of extr do not have a precise
			% meaning, but they are simple scale factors determinimi the
			% extension of the null space to be considered for the
			% discretization (they stretch'' the unit vector defining the
			% direction of the null space, G_red).

			extr = [NaN, NaN];
	end
	% Number of points in the discretization of the null space
	N_NS = 101;
	% Discretize null space
	NS_red = t_red_DI + G_red'.*repmat(linspace(extr(1), extr(2), N_NS)', 1, A);
	% Null space uncertainty, method by Palací-López et al. (2019), Chemometr. Intell. Lab. Syst. 194, 103848
	NS_red_cl = pls_inversion_null_space_uncertainty_nonlin(NS_red, Q_tilde_red, y_targ_red_IbO, alpha, sigma_sq, SSE, dof, N, V_Y - sum(y_free));
	% Confidence limits of the null space
	NS_red_lcl = NS_red - NS_red_cl;
	NS_red_ucl = NS_red + NS_red_cl;
else
	% No null space
	G_red = [];
	t_red_DI_lcl = [];
	t_red_DI_ucl = [];
	NS_red = [];
	NS_red_lcl = [];
	NS_red_ucl = [];
end

%% Plots

%---> REVIEW OR MODIFY SETTINGS FOR PLOTS <---%
% Note that plots are generated only for the example datasets provided with this
% script. If custom data are used, no plot is produced.

% Colour sequence
colours = [
	0.3000,	0.3000,	0.3000	% Data
	0.8392, 0.1529, 0.1569	% Target (red)
	0.1216, 0.4667, 0.7059	% DI (blue)
	1.0000, 0.4980, 0.0549	% RDI (orange)
	0.1725, 0.6275, 0.1725	% IbO (green)
	0.0000, 0.0000, 0.0000	% Null space (balck)
];

% Marker styles
lsm = {
	'x'	% Data
	'h'	% Target
	'o'	% DI
	's'	% RDI
	'^'	% IbO
};

% Maker size
ms = 6;
% Marker size scale factors
mssf = [
	1		% Data
	1.5		% Target
	1		% Solutions
];

% Check if plots are o be produced
if strcmp(datasel, 'palaci-lopez') || strcmp(datasel, 'bioreactor')
	
	% Preliminary operations
	switch datasel
		case 'palaci-lopez'
			figure(1)
			clf
			tiledlayout(2, 10, 'Padding', 'none', 'TileSpacing', 'tight')
		case 'bioreactor'
			figure(1)
			clf
			tiledlayout(2, 8, 'Padding', 'none', 'TileSpacing', 'tight')
	end
	
	% Quality space
	switch datasel
		case 'palaci-lopez'
			figure(1)
			nexttile(1, [1, 3])
			hold on
			plot(1:N, Y(:, 1), lsm{1}, 'MarkerSize', ms*mssf(1), 'MarkerEdgeColor', colours(1, :))
			plot(N + 1, y_des(:, 1), lsm{2}, 'MarkerSize', ms*mssf(2), 'MarkerEdgeColor', colours(2, :), 'MarkerFaceColor', colours(2, :)*0.7 + [1, 1, 1]*0.3)
			plot(N + 2, y_des_DI(:, 1), lsm{3}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(3, :), 'MarkerFaceColor', colours(3, :)*0.7 + [1, 1, 1]*0.3)
			plot(N + 3, y_des_RDI(:, 1), lsm{4}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(4, :), 'MarkerFaceColor', colours(4, :)*0.7 + [1, 1, 1]*0.3)
			plot(N + 4, y_des_IbO(:, 1), lsm{5}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(5, :), 'MarkerFaceColor', colours(5, :)*0.7 + [1, 1, 1]*0.3)
			hold off
			title('Quality space')
			xlabel('Observation')
			ylabel('Y_1')
		case 'bioreactor'
			figure(1)
			nexttile(1, [1, 2])
			hold on
			plot(Y(:, 1), Y(:, 2), lsm{1}, 'MarkerSize', ms*mssf(1), 'MarkerEdgeColor', colours(1, :))
			plot(y_des(:, 1), y_des(:, 2), lsm{2}, 'MarkerSize', ms*mssf(2), 'MarkerEdgeColor', colours(2, :), 'MarkerFaceColor', colours(2, :)*0.7 + [1, 1, 1]*0.3)
			plot(y_des_DI(:, 1), y_des_DI(:, 2), lsm{3}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(3, :), 'MarkerFaceColor', colours(3, :)*0.7 + [1, 1, 1]*0.3)
			plot(y_des_RDI(:, 1), y_des_RDI(:, 2), lsm{4}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(4, :), 'MarkerFaceColor', colours(4, :)*0.7 + [1, 1, 1]*0.3)
			plot(y_des_IbO(:, 1), y_des_IbO(:, 2), lsm{5}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(5, :), 'MarkerFaceColor', colours(5, :)*0.7 + [1, 1, 1]*0.3)
			hold off
			title('Quality space')
			xlabel('Y_1')
			ylabel('Y_2')
	end
	
	% Designed process conditions
	switch datasel
		case 'palaci-lopez'
			figure(1)
			tilemap = [11, 12, 13, 14, 15];
			for v = 1:V_X
				nexttile(tilemap(v))
				hold on
				plot(0, X(:, v), lsm{1}, 'MarkerSize', ms*mssf(1), 'MarkerEdgeColor', colours(1, :))
				plot(0, x_des(:, v), lsm{2}, 'MarkerSize', ms*mssf(2), 'MarkerEdgeColor', colours(2, :), 'MarkerFaceColor', colours(2, :)*0.7 + [1, 1, 1]*0.3)
				plot(0, x_des_DI(:, v), lsm{3}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(3, :), 'MarkerFaceColor', colours(3, :)*0.7 + [1, 1, 1]*0.3)
				plot(0, x_des_RDI(:, v), lsm{4}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(4, :), 'MarkerFaceColor', colours(4, :)*0.7 + [1, 1, 1]*0.3)
				plot(0, x_des_IbO(:, v), lsm{5}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(5, :), 'MarkerFaceColor', colours(5, :)*0.7 + [1, 1, 1]*0.3)
				set(gca, 'XTick', [])
				ylabel(sprintf('X_%d', v))
			end
		case 'bioreactor'
			figure(1)
			tilemap = [9, 10, 11, 12];
			for v = 1:V_X
				nexttile(tilemap(v))
				hold on
				plot(0, X(:, v), lsm{1}, 'MarkerSize', ms*mssf(1), 'MarkerEdgeColor', colours(1, :))
				plot(0, x_des(:, v), lsm{2}, 'MarkerSize', ms*mssf(2), 'MarkerEdgeColor', colours(2, :), 'MarkerFaceColor', colours(2, :)*0.7 + [1, 1, 1]*0.3)
				plot(0, x_des_DI(:, v), lsm{3}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(3, :), 'MarkerFaceColor', colours(3, :)*0.7 + [1, 1, 1]*0.3)
				plot(0, x_des_RDI(:, v), lsm{4}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(4, :), 'MarkerFaceColor', colours(4, :)*0.7 + [1, 1, 1]*0.3)
				plot(0, x_des_IbO(:, v), lsm{5}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(5, :), 'MarkerFaceColor', colours(5, :)*0.7 + [1, 1, 1]*0.3)
				set(gca, 'XTick', [])
				ylabel(sprintf('X_%d', v))
			end
	end
	
	% Score space
	switch datasel
		case 'palaci-lopez'
			figure(1)
			nexttile(6, [2, 5])
		case 'bioreactor'
			figure(1)
			nexttile(5, [2, 4])
	end
	hold on
	a = 1;
	if ~isempty(G_red)
		hns = plot(t_red_DI(a) + G_red(a, :).*extr, t_red_DI(a + 1) + G_red(a + 1, :).*extr, '-', 'Color', colours(6, :));
		hnslcl = plot(t_red_DI_lcl(a) + G_red(a, :).*extr, t_red_DI_lcl(a + 1) + G_red(a + 1, :).*extr, '--', 'Color', colours(6, :));
		hnsucl = plot(t_red_DI_ucl(a) + G_red(a, :).*extr, t_red_DI_ucl(a + 1) + G_red(a + 1, :).*extr, ':', 'Color', colours(6, :));
		plot(NS_red_lcl(:, a), NS_red_lcl(:, a + 1), '--', 'Color', colours(6, :))
		plot(NS_red_ucl(:, a), NS_red_ucl(:, a + 1), ':', 'Color', colours(6, :))
	else
		hns = plot(NaN, NaN, '-', 'Color', colours(6, :));
		hnslcl = plot(NaN, NaN, '--', 'Color', colours(6, :));
		hnsucl = plot(NaN, NaN, ':', 'Color', colours(6, :));
	end
	h0 = plot(T(:, a), T(:, a + 1), lsm{1}, 'MarkerSize', ms*mssf(1), 'MarkerEdgeColor', colours(1, :));
	h1 = plot(t_des(:, a), t_des(:, a + 1), lsm{2}, 'MarkerSize', ms*mssf(2), 'MarkerEdgeColor', colours(2, :), 'MarkerFaceColor', colours(2, :)*0.7 + [1, 1, 1]*0.3);
	h2 = plot(t_des_DI(:, a), t_des_DI(:, a + 1), lsm{3}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(3, :), 'MarkerFaceColor', colours(3, :)*0.7 + [1, 1, 1]*0.3);
	h3 = plot(t_des_RDI(:, a), t_des_RDI(:, a + 1), lsm{4}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(4, :), 'MarkerFaceColor', colours(4, :)*0.7 + [1, 1, 1]*0.3);
	h4 = plot(t_des_IbO(:, a), t_des_IbO(:, a + 1), lsm{5}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(5, :), 'MarkerFaceColor', colours(5, :)*0.7 + [1, 1, 1]*0.3);
	[x_cl, y_cl] = ellipse(0, 0, l_T(a), l_T(a + 1));
	plot(x_cl, y_cl, 'r-.')
	axis equal
	xlims = xlim;
	ylims = ylim;
	plot(xlims, [0, 0], '-', 'Color', [0.8, 0.8, 0.8])
	plot([0, 0], ylims, '-', 'Color', [0.8, 0.8, 0.8])
	handles = findobj(gca, 'Type', 'Line');
	uistack(handles(1:3), 'bottom');
	hold off
	xlim(xlims);
	ylim(ylims);
	title(sprintf('Scores plot: LVs %d and %d%', a, a + 1))
	xlabel(sprintf('T_{%d} (EV_X = %.4g%%, PV_Y = %.4g%%)', a, EV_X(a)*100, PV_Y(a)*100))
	ylabel(sprintf('T_{%d} (EV_X = %.4g%%, PV_Y = %.4g%%)', a + 1, EV_X(a + 1)*100, PV_Y(a + 1)*100))
	legend([h0, h1, h2, h3, h4, hns, hnslcl, hnsucl], {'Data', 'Target', 'DI', 'RDI', 'IbO' 'Null space', 'NS LCL', 'NS UCL'}, 'Location', 'northeast')
	
	% Diagnostic plot
	switch datasel
		case 'palaci-lopez'
			figure(1)
			nexttile(4, [1, 2])
		case 'bioreactor'
			figure(1)
			nexttile(3, [1, 2])
	end
	hold on
	a = 1;
	plot(T_sq, SRE, lsm{1}, 'MarkerSize', ms*mssf(1), 'MarkerEdgeColor', colours(1, :))
	plot(T_sq_DI, SRE_DI, lsm{3}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(3, :), 'MarkerFaceColor', colours(3, :)*0.7 + [1, 1, 1]*0.3)
	plot(T_sq_RDI, SRE_RDI, lsm{4}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(4, :), 'MarkerFaceColor', colours(4, :)*0.7 + [1, 1, 1]*0.3)
	plot(T_sq_IbO, SRE_IbO, lsm{5}, 'MarkerSize', ms*mssf(3), 'MarkerEdgeColor', colours(5, :), 'MarkerFaceColor', colours(5, :)*0.7 + [1, 1, 1]*0.3)
	tempa = plot(xlim, lim_SRE*[1, 1], 'r-.');
	tempb = plot(lim_T_sq*[1, 1], ylim, 'r-.');
	xlims = xlim;
	ylims = ylim;
	delete(tempa)
	delete(tempb)
	plot(xlims, lim_SRE*[1, 1], 'r-.')
	plot(lim_T_sq*[1, 1], ylims, 'r-.')
	hold off
	xlim(xlims);
	ylim(ylims);
	title(sprintf('SRE vs. T^2 with %d LVs', A))
	xlabel('T^2')
	ylabel('SRE')
end
