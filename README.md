# Streaming, Memory-Limited, Federated PCA Revisited!

In this work, we present a novel federated algorithm for PCA that 
is able to adaptively estimate the rank `r` of the dataset and compute 
its r-leading principal components when only finite memory, 
specifically `O(dr)`, is available. This inherent adaptability implies 
that the rank `r` does not have to be supplied as a fixed hyper-parameter 
which is beneficial when the underlying data distribution is not known 
in advance - such as in a streaming setting. Numerical simulations show 
that, while using limited-memory, our algorithm exhibits state-of-the-art 
performance that closely matches or outperforms traditional non-federated 
algorithms, and in the absence of communication latency, it exhibits 
attractive horizontal scalability.

# Requirements

The code is generally self contained and all datasets are included or 
generated thus, in theory, just having `Matlab`  installed should be more 
than enough. It has to be noted though that due the recent `Matlab` changes 
on how it handles character and string arrays you should use a recent 
version of it -- the code was developed and tested in `Matlab` `2019a` 
build `9.6.0.1099231` but was tested also on versions `2018a` and `2018b`; 
moreover, to address different OSes, care has been taken so that this 
code runs without any problems both on Windows-based machines as well 
as Unix-based ones.

# Comparisons

In this instance we perform comparisons using both synthetic and real 
datasets against a few similar methods which compute in part or fully an 
approximate *memory-limited, streaming r-truncated PCA*. To make the 
comparison fair we perform single node experiments for each method and
compare their respective outputs. 

 * Federated PCA (https://arxiv.org/abs/1907.08059)
 * Power Method (https://arxiv.org/pdf/1307.0032.pdf)
 * Frequent Directions (https://arxiv.org/abs/1501.01711.pdf)
 * Robust Frequent Directions (https://arxiv.org/pdf/1705.05067.pdf)
 * GROUSE (https://arxiv.org/pdf/1702.01005.pdf)
 * SPIRIT (https://dl.acm.org/citation.cfm?id=1083674)

Note that the rank adjusting experiments are performed using 
only SPIRIT against our method as is the only method that 
has an explicit rank estimation mechanism via energy 
thresholding. Finally, note that S(A)PCA is in spirit 
similar to [MOSES][5] and inherits most of its properties 
and thus no comparison is made against it.

# Running the comparison

Running the comparison is simple -- just `cd` to the cloned 
`federated_pca`  directory within `Matlab` and then run the 
respective test files - brief explanation of what they do
is shown below:

 * [`test_sapca_real.m`](test_sapca_real.m): tests for the real datasets.
 * [`test_sapca_synthetic.m`](test_sapca_synthetic.m): tests for the synthetic datasets.
 * [`test_sapca_federated.m`](test_sapca_federated.m): performs the federated tests.
 * [`test_subspace_merge_error.m`](test_subspace_merge_error.m): performs the subspace merging error tests.
 * [`test_time_order.m`](test_time_order.m): performs the time order invariance tests.

Please note that you can tweak the relevant section values 
if you want to run slightly different experiments but if 
you want to reproduce the results in the paper please leave 
these values as-is.

# Synthetic Datasets

The synthetic dataset is measured using random vectors drawn from a power
law distribution with the following alpha values in this instance: `0.0001`, `0.001`, 
`0.5`, `1`, `2` and `3` while lambda always set to `1`. Practically speaking
this is eloquently materialised by using the following segment:

```Matlab
% generate the singular spectrum
dd = lambda*(1:n).^(-alpha);
% generate Sigma
Sigma = diag(dd);
% random initialization of S basis
S = orth(randn(n));
% given S and Sigma generate the dataset (Y)
Y = (S * Sigma * randn(n, T))/sqrt(T-1);
```

# Real Datasets

The real datasets are the the ones supplied with [this][3] paper 
retrieved from [here][1] and they are the following:

 * Light Data (48x7712)
 * Humidity Data (48x7712)
 * Volt Data (46x7712)
 * Temperature Data (56x7712)

# Error metrics

To compare S(A)PCA against Power Method, FD/RFD, and GROUSE we employ the 
following two metrics:

 * The Frobenius norm of `Yr` vs `Y` columns seen so far normalised using 
    their respective arrival time.
 * The final Frobenius norm of `Yr` vs `Y` normalised by the final `T`.

## Normalised Frobenius norm over time normalised with current T

The error metrics are calculated using the Frobenius norm for the
matrix columns seen so far normalised by the current time. The full 
formula to find the error at column `k` would be:

```Matlab
ErrFro(k) = sum(sum((Y(:, 1:t)-YrHat_c).^2, 1))/t;
```

Where `YrHat_c` is:

```Matlab
SrHatTemp = SrHat(:, 1:r); % r-truncation of the SVD
% SrHat in this instance is the previous block subspace estimation
YrHat_c = (SrHat*SrHat')*Y(:, 1:k*B); 
```

## MSE of the final Subspace vs Offline  

The other metric is the MSE between the subspace produced by an 
offline `PCA` and the one approximated by each method. This enables
us to see how each approximation differs from the target objective.

## Principal Components vs Error over time

For the adaptive methods (ours and SPIRIT) we also test how the 
evolution of Principal Components (PC's) occurs over time with 
respect to the errors previously mentioned - ideally, we'd like 
to have the lowest error possible with the fewest Principal 
Components.


# Federated Tests

In order to check how the algorithm would perform in a federated 
setting we construct a tree hierarchy which comprises out of
aggregators and edge nodes which are responsible for merging
and PCA computation respectively. We report the actual
and amortised execution speed results for various depths 
and dataset sizes.

# Plots

A number of plots are generated while running the comparison and for
convenience they are printed into a generated directory under the 
`graph` directory. Each directory is named using the current 
timestamp upon creation as its name and the timestamp format 
follows the [ISO-8601][4] standard.

Additionally, the printing function is flexible enough to able to export 
in three commonly used formats concurrently -- namely `png`, `pdf`, and 
`fig` for easier processing. Of course, by toggling the appropriate flags 
printing to `pdf` and `fig` can be disabled thus saving space. For 
brevity these are the following:

```MatLab
% printing flags
pflag = 1;              % print resulting figures to ./graphs/
pdf_print = 0;          % print resulting figures as .pdf
fig_print = 1;          % print resulting figures as .fig
```

Please note that `Matlab` is sometimes picky when exporting `pdf` 
figures on high-dpi displays... so your mileage may vary!

# Code Organisation

The code is self-contained and a brief explanation of what each file does follows. The 
files are ordered in (descending) lexicographical order:

 * `fd.m`: Implementation of Frequent Directions.
 * `fd_rotate_sketch.m`: helper method for both Frequent Directions methods.
 * `fdr.m`: Implementation of Robust Frequent Directions.
 * `grams.m`: Gram-Schmidt orthogonalization for a given matrix.
 * `grouse.m`: Original `GROUSE` algorithm code as provided from its authors.
 * `merge_subspaces.m`: Merge two subspaces using different techniques.
 * `mitliag_pm.m`: Implementation of Mitliagkas Power Method for Streaming PCA.
 * `my_grouse.m`: Wrapper to run `grouse.m` which sets execution parameters (as seen [here][2]).
 * `my_toc.m`: function that processes the `toc` with better formatting.
 * `print_fig.m`: Prints figures in different formats (i.e.: `pdf`, `png`, and `fig`).
 * `README.md`: This file, a brief README file.
 * `real_sapca_eval.m`: runs the evaluation for a provided (real) dataset.
 * `sapca_edge`: directly and incrementally computes SA-PCA within each edge node.
 * `setup_vars.m`: sets up the environment variables.
 * `spca_edge.m`: directly and incrementally computes S-PCA within each edge node.
 * `spectrum_adaptive.m`: plot helper for the singular value approximation vs ground truth.
 * `SPIRIT.m`: Original `SPIRIT` algorithm as provided from its authors (as seen [here][1]).
 * `synthetic_data_gen.m`: function which generates a matrix with random vectors from a power law distribution.
 * `test_sapca_real.m`: performs the tests for the real datasets - which are provided.
 * `test_sapca_synthetic.m`: performs the tests for the synthetic datasets - which are generated.
 * `test_sapca_federated.m`: performs the federated tests.
 * `test_subspace_merge_error.m`: performs the subspace merging error tests and was used as a test-bed.
 * `test_time_order.m`: performs the time order invariance tests.
 * `updateW.m`: helper function for SPIRIT, performs the update of the subspace for each datapoint.


# License

This code is licensed under the terms and conditions of GPLv3 unless otherwise stated. 
The actual paper is governed by a separate license and the paper authors retain their 
respective copyrights.

# Acknowledgement

If you find our paper useful or use this code, please consider citing our work as such:

```
@misc{1907.08059,
Author = {Andreas Grammenos and Rodrigo Mendoza-Smith and Cecilia Mascolo and Jon Crowcroft},
Title = {Federated PCA with Adaptive Rank Estimation},
Year = {2019},
Eprint = {arXiv:1907.08059},
}
```

# Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

[1]: http://www.cs.cmu.edu/afs/cs/project/spirit-1/www/
[2]: http://web.eecs.umich.edu/~girasole/grouse/
[3]: http://www.cs.albany.edu/~jhh/courses/readings/desphande.vldb04.model.pdf
[4]: https://en.wikipedia.org/wiki/ISO_8601
[5]: https://github.com/andylamp/moses