# Description of functions:

"Jisst_single": Single-factor JisstPCA. The signal in $\mathcal{X}, \mathcal{Y}$ are scaler $d_{x}$ and $d_{y}$.

"Jisst_multi": Multi-factor JisstPCA. The signals in each layer of tensor are scalers.

"dJisst_single": Single-factor JisstPCA with diagonal matrix $D_{x}, D_{y}$ as signals.

"dJisst_multi": Multi-factor JisstPCA. The signals in each layer of tensor are diagonal matrices $D_{x}, D_{y}$. 

"iHOSVD": Integrated high-order singular value decomposition.

"iHOOI": Integrated high-order orthogonal iteration.

"bic_sst": Value of BIC, given $\mathcal{X}, \mathcal{Y}$ and input ranks $r_{x}, r_{y}$. This is a function of $(r_{x}, r_{y})$.

"bic_def": Rank selection based on finding $(r_{x}, r_{y})$ that minimize the BIC value, for input data $\mathcal{X}$ and $\mathcal{Y}$.

"init": Return the spectral initialization based on $\mathcal{X}, \mathcal{Y}$.

And other functions used in the main functions are:

"tr_prod": Trace product between a semi-symmetric tensor and a low-rank matrix.

"gtr_prod": Generalized trace product between a semi-symmetric tensor, diagonal matrix and a low-rank matrix.

