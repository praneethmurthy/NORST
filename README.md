This folder contains the code accompanying paper.

	[1] "Nearly Optimal Robust Subspace Tracking", Praneeth Narayanamurthy and Namrata Vaswani, ICML, 2018 (long version available at https://arxiv.org/abs/1712.06061).

List of main files:

	1. DemoFB.m - Wrapper containing the real video foreground background separation. 
	2. DemoDynRPCA.m - Wrapper containing the simulated data experiments. 
	3. NORST - main function which implements the NORST algorithm for subspace tracking.
	4. Offline_NORST - main function which implements the Offline Norst variant for dynamic (and static) Robust PCA
	5. NORST_real - main function which implements foreground-background separation.
	6. phase_transition_generation.m - file to generate phase transition plots
	7. xmin_variation.m - file to verify the effect of varying x_min. 

Folders:

	1. YALL1 - folder containing files to implement ell-1 minimization.
	2. PROPACK - Linear Algebra toolbox for MATLAB

Other files:

	1. ncrpca -- code implemented Non-convex Robust PCA, NIPS 14 downloaded from authors' website and its accompaniments lansvd, lanbpro etc
	2. cgls -- fast method to implement least squares


For any further questions/suggestions please contact me @ pkurpadn iastate edu (insert obvious symbols)

