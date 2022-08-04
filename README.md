# GRAND-MATLAB

Ken R. Duffy, 2018-2022.

Guessing Random Additive Noise Decoding (GRAND)
For further details see: https://www.granddecoder.mit.edu/

Non-parallelized MATLAB implementations of: GRAND (hard detection); basic ORBGRAND (soft detection); 1-line ORBGRAND (soft detection).

Subject to license:
GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

For an [n,k] code, where k information bits become n coded bits,
GRAND algorithms accurately and efficiently decode codes where n-k
is moderate. This MATLAB implementation is solely intended to be
instructive and highly-parallelised implementations are possible.
As a result, obtaining the full performance a code with n-k>16 will
be time-consuming with the present implementation.

The following should be cited in association with results from this code.

K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random 
additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 
4023–4040, 2019.

K. R. Duffy, “Ordered reliability bits guessing random additive noise 
decoding," in IEEE ICASSP, 2021, pp. 8268–8272. 

K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
random additive noise decoding,” arXiv:2202.13951, 2022.
