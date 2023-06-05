# Description of functions:

## Main functions of JisstPCA:

"JisstPCA": Main function of joint-integrated semi-symmetric tensor PCA algorithm. 

<ul>
  <li>Required input:</li>
  <ul>
    <li> X, Y: Two semi-symmetric tensors. </li>
    <li> K: Number of layers that users care. If $K = 1$ then the function degenerates to single-factor JisstPCA. </li>
  </ul>
  
  <li> Optional inputs (if users do not specify them then JisstPCA will apply the default values): </li>
  <ul>
    <li> u0: Initialization. (Default: spectral initialization.) </li>
    <li> rx, ry: Ranks of $\mathcal{X}$ and $\mathcal{Y}$. (Default: BIC + deflation for rank selection.) </li>
    <li> lambda: Relative importance of $\mathcal{X}$ compared with $\mathcal{Y}$, which is a vector of length $K$. (Default: $\frac{\|\mathcal{X}\|}{\|\mathcal{X}\|+\|\mathcal{Y}\|}$ for each element.) </li>
    <li> tol: Tolerance level. (Default: tol = 0.0001.) </li>
    <li> max_iter: Maximum iteration number. (Default: max_iter = 100.) </li>
    <li> deflation: Deflation strategy: subtract deflation, project deflation or partial project deflation (see JisstPCA.m for more detail). (Default: deflation = 0, subtract deflation.) </li>
    <li> rank_max: When using BIC for rank selection, the largest possible rank. (Default: rank_max = 5.) </li>
    <li> method: When using BIC for rank selection, method = 1 assumes $\mathcal{X}$ and $\mathcal{Y}$ are of differnet ranks, while method = 2 assumes $\mathcal{X}$ and $\mathcal{Y}$ are of the same ranks. (Default: method = 1.) </li>
  </ul>
  
  <li> Outputs: </li>
  <ul>
    <li> u_est, V_est, W_est: Cells with $K$ elements, with $k$-th element being the estimation of $\boldsymbol{u}_{k}, \boldsymbol{V}_{k}$ or $\boldsymbol{W}_{k}$. </li>
    <li> d_est: Matrix in $\mathbb{R}^{2 \times K}$. The first row is estimation of signal of $\mathcal{X}$, while the second row is estimation of signal of $\mathcal{Y}$. </li>
  </ul>
  
  <li> Note and examples: </li>
  <ul>
    <li> Note: When using this function, users have to specify the name if they have their own choice of optional inputs! </li>
    <li> Example 1: If users only have $\mathcal{X}, \mathcal{Y}$ and $K$, then they can implement JisstPCA as: result_1 = JisstPCA(X, Y, K). </li>
    <li> Example 2: If users also have $r_{x}$ and $r_{y}$, then this function should be used as: result_2 = JisstPCA(X, Y, K, 'rx', rx, 'ry', ry). </li>
  </ul>
</ul>

"dJisstPCA": Main function of joint-integrated semi-symmetric tensor PCA algorithm with signal being diagonal matrix. The usage is similar as "JisstPCA" and the only difference is the output of "dJisstPCA" includes estimation of diagonal matrices signals.

More example and detail interpretion of JisstPCA/dJisstPCA can be found in "example" folder and JisstPCA.m/dJisstPCA files.

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

