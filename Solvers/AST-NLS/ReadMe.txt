************************************************************************
            AST-NLS: Image Denoising via Adaptive Soft-Thresholding Based on Non-Local Samples (beta version)
************************************************************************


Usage
====================

AST-NLS is called in the following way:

               	% clean_image = double( imread( 'image.png' ) ); 
		% noisy_image = clean_image + noise_sigma * randn( size ( clean_image ) );
		  denoised_image = ast_nls( noisy_image, noise_sigma );


Requirements
====================

		MS Windows;
		Matlab 2010b or later versions.
