# this file contains collections of proxes we learned in the class
import numpy as np
from scipy.optimize import bisect

# =============================================================================
# TODO Complete the following prox for simplex
# =============================================================================

# Prox of capped simplex
# -----------------------------------------------------------------------------
def prox_csimplex(z, k):
    """
    Prox of capped simplex
	    argmin_x 1/2||x - z||^2 s.t. x in k-capped-simplex.

	input
	-----
    z : arraylike
		reference point
	k : int
		positive number between 0 and z.size, denote simplex cap

	output
	------
    x : arraylike
		projection of z onto the k-capped simplex
        """
    # safe guard for k
    assert 0<=k<=z.size, 'k: k must be between 0 and dimension of the input.'
	# TODO do the computation here
	# Hint: 1. construct the scalar dual object and use `bisect` to solve it.
	#		2. obtain primal variable from optimal dual solution and return it.
	#
    
    one = np.ones(z.size)
    dual = lambda y: np.sum(np.maximum(np.minimum(z-y*one,one),np.zeros(z.size)))-k
    y = bisect(dual,np.min(z)-1,np.max(z))
    x = np.maximum(np.minimum(z-y*one,one),np.zeros(z.size))
    return x

