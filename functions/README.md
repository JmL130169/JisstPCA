# Description of functions:

## Main functions of JisstPCA:

## Other main functions: 

"Jisst_single": Single-factor JisstPCA. The signal in $\mathcal{X}, \mathcal{Y}$ are scaler $d_{x}$ and $d_{y}$.

"dJisst_single": Single-factor JisstPCA with diagonal matrix $D_{x}, D_{y}$ as signals.

"iHOSVD": Integrated high-order singular value decomposition.

"iHOOI": Integrated high-order orthogonal iteration.

"bic_sst": Value of BIC, given $\mathcal{X}, \mathcal{Y}$ and input ranks $r_{x}, r_{y}$. This is a function of $(r_{x}, r_{y})$.

"bic_dsst": The counterpart of "bic_sst" for diagonal signal case.

"bic_def": Rank selection based on finding $(r_{x}, r_{y})$ that minimize the BIC value, for input data $\mathcal{X}$ and $\mathcal{Y}$.

"bic_diag": The counterpart of "bic_def" for diagonal signal case.

"bic_def_1" ("bic_def_2"): The output of estimating $r_{x}$ ($r_{y}$) in "bic_def" (this function is used in main function "JisstPCA").

"bic_diag_1" ("bic_diag_2"): The output of estimating $r_{x}$ ($r_{y}$) in "bic_diag" (this function is used in main function "dJisstPCA").

"init": Return the spectral initialization based on $\mathcal{X}, \mathcal{Y}$.

## Other functions used in the main functions are:

"tr_prod": Trace product between a semi-symmetric tensor and a low-rank matrix.

"gtr_prod": Generalized trace product between a semi-symmetric tensor, diagonal matrix and a low-rank matrix.

"sin_do": $\sin\Theta$ distance between two subspaces spanned by two matrices.

