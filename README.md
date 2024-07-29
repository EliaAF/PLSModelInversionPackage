# PLS Model Inversion Package
#### (C) Elia Arnese Feffin – February 20, 2024

Version: 1.0.0

Date: 2024-02-20

Author: Elia Arnese Feffin elia.arnesefeffin@phd.unipd.it/elia249@mit.edu

## Package description, requisites, and usage

This repository contains MATLAB code for calibration, selection, application, and inversion of partial least-squares (PLS) regression models. The methods are implemented as a collection of MATLAB functions. All functions have been developed and tested on MatLab R2022a for MacOS and are meant to be used with a minimal MATLAB installation. The required MATLAB toolboxes are the following.

* Statistics and Machine Learning Toolbox
* Optimization Toolbox

All functions are extensively documented in order to clarify their rationale, what inputs and outputs are, and what options are available. Options are to be passed with the key-value logic. The documentation is accessible using the help command in MATLAB. As functions have a high degree of interdependency, it is recommended to keep all of them in the current working directory.

The file `LVMI_case_studies.m` serves as an example on how to use the main functions of this package. The file is a MATLAB script guiding the user through the workflow for PLS model calibration, selection, and inversion. Three methods for model inversion (direct inversion, regularized direct inversion, and inversion by optimization) are carried out and compared. Inputs to be reviewed and/or set by the user are distributed throughout this script. As inputs are not gathered in a single section of the script and are spread throughout the code, it is warmly recommended to execute the code section by section, carefully reviewing all the instructions in doing so. Detailed instructions on how to use the script are provided in the “Instructions and information” section of the script itself.

The script offers two pre-implemented case studies.

* The “example data” provided by Palací-López et al. (2019), Chemometr. Intell. Lab. Syst. 194, 103848.
* The “bioreactor case study” discussed by Arnese-Feffin et al. (2022), Chemometr. Intell. Lab. Syst. 227, 104610.

Additionally, the script allows the user to provide custom data for PLS model calibration and inversion. In this case, please carefully review the instructions spread throughout the script. Bear in mind that this script is to be interpreted as an example rather than a general interface.

## Package content

 A list of the functions provided, with a brief description of their rationales, is reported below.
 
* Primary functions
	- `LVMI_case_studies.m`: Primary script for PLS model calibration and inversion
	- `build_pls.m`: Primary function for PLS model calibration
	- `apply_pls.m`: Primary function for PLS model application
	- `cross_validate_pls.m`: Primary function for PLS model cross-validation
	- `cross_validate_pls_nested.m`: Primary function for PLS model nested
	  cross-validation
* Fundamental PLS modeling
	- `pls_calibration_nipals.m`: PLS model calibration by the nonlinear iterative PLS algorithm
	- `pls_calibration_simpls.m`: PLS model calibration by the statistically inspired modification of the PLS algorithm
	- `pls_prediction.m`: Simple application of calibrated PLS models
	- `pls_prediction_uncertainty.m`: Uncertainty in prediction of PLS models
	- `pls_dof.m`: degrees of freedom of the PLS model
	- `initialize_pls.m`: Initialization of the PLS model structure
* PLS model diagnostics
	- `T_sq_statistic.m`: Hotelling's T<sup>2</sup> statistic
	- `Q_statistic.m`: Q statistic, also called squared-reconstruction (or prediction) error statistic
	- `F_limit.m`: Control limit by the F distribution approach
	- `jackson_mudholkar_limit.m`: Control limit by the Jackson & Mudholkar approach
	- `chi_sq_limit.m`: Control limit by the χ<sup>2</sup> distribution approach
	- `kde_limit.m`: Control limit by the kernel density estimation approach
* PLS model inversion
	- `pls_inversion_DI.m`: PLS model direct inversion
	- `pls_inversion_RDI.m`: PLS model regularized direct inversion
	- `pls_inversion_IbO.m`: PLS model inversion by optimization
	- `pls_inversion_null_space_uncertainty_lin.m`: Null space uncertainty with the method by Facco et al. (2015), Ind. Eng. Chem. Res. 54, 5128-5138
	- `pls_inversion_null_space_uncertainty_nonlin.m`: Null space uncertainty with the method by Palací-López et al. (2019), Chemometr. Intell. Lab. Syst. 194, 103848
	- `pls_IBO_objective.m`: Objective function for PLS IbO
	- `pls_IBO_initialize_constraints.m`: Initialization of the structure for
	  linear constraints in PLS IbO
	- `pls_IBO_linear_constraints.m`: Assemble linear constraint matrices for PLS IbO
	- `pls_IBO_nonlinear_constraints_SD.m`: Nonlinear constraint function on squared distances (T<sup>2</sup> or Q) in PLS IbO
	- `pls_IBO_nonlinear_constraints_NS.m`: Nonlinear constraint function on the null space in PLS IbO
	- `pls_IBO_nonlinear_constraints_wrapper.m`: Wrapper function for nonlinear constraints in PLS IbO
* Utilities
	- `autoscale.m`: Autoscaling of data matrices
	- `scale_by.m`: Scaling of a data matrix with given means and standard		  deviations
	- `rescale_by.m`: Rescaling of a data matrix to given means and standard		  deviations
	- `slim_cv_pls.m`: Efficient PLS cross-validation with the leave-one-out observation splitting scheme
	- `cross-validation_grouping.m`: Splitting schemes for observations in cross-validation
	- `information_criteria.m`: Collection of information criteria for model selection
	- `prob_dist_IC.m`: Probability distributions of the information criteria
	- `kernel_density_estimation.m:` Efficient univariate kernel density estimation
	- `gaussian_kernel.m`: Gaussian kernel function for univariate kernel density estimation
	- `ellipse.m`: Plot of an ellipse given the center and semiaxes
	- `reg_inv.m`: Regularized matrix inversion
* Datasets
	- `dataset_palacilopez2019_ex.mat`: Example data from Palací-López et al. (2019), Chemometr. Intell. Lab. Syst. 194, 103848
	- `dataset_batchcellculture.mat`: Bioreactor case-study data by Arnese-Feffin et al. (2022), Chemometr. Intell. Lab. Syst. 227, 104610
* Other files
	- `README.md`: This file
	- `GPL-3.0.txt`: GNU General Public License version 3

## Attribution

To attribute credit to the author of the software, please refer to the following Journal Paper.

E. Arnese-Feffin, P. Facco, F. Bezzo, and  M. Barolo (2024): Systematizing product design by latent-variable modeling ‒ A unifying framework for the formulation and solution of PLS model-inversion problems. _Chemical Engineering Science_ **299**, 120505. DOI: [10.1016/j.ces.2024.120505](https://doi.org/10.1016/j.ces.2024.120505).

## License agreement

All the files provided are covered by the GNU General Public License version 3 (GPL-3.0); a copy of the licence can be found in this folder (GPL-3.0.txt), or online at https://www.gnu.org/licenses/gpl-3.0.html. This licence protects the code as open-source. Key points of GPL-3.0 are as follows.

* Attribution of credit to the author is required if the software is used.
* Free use of the software is allowed, even for commercial purposes.
* Redistribution of the software for commercial purposes is prevented, as any redistribution must be released under the GPL-3.0 licence, therefore as a free and open-source software.
