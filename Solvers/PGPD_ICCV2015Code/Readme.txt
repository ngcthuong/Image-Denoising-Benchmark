% ===============================================================

The code in this package implements the Patch Group Prior Denoising (PGPD) method for
image denoising as described in the following paper:

  

Jun Xu, Lei Zhang, Wangmeng Zuo, David Zhang, and Xiangchu Feng
Patch Group Based Nonlocal Self-Similarity Prior Learning for Image Denoising
 
IEEE Int. Conf. Computer Vision (ICCV), Santiago, Chile, December 2015.



Please cite the paper if you are using this code in your research.

Please see the file License.txt for the license governing this code.

 
 
Version:       
1.0 (04/09/2015), see ChangeLog.txt
  

Contact:       
Jun Xu <csjunxu@comp.polyu.edu.hk>

% ===============================================================


Overview

------------

The code for learning Patch Group Prior is implemented in the folder "PG-GMM_Traning", which relies
on the training images in the subfolder "Kodak24".

The function "Demo_denoising" demonstrates denoising with the learned Patch Group Prior models 
introduced in the paper, which can all be found in the folder "model".



Dependency

------------

This code is implemented purely in Matlab2014b and doesn't depends on any other toolbox.



Contact

------------

If you have questions, problems with the code, or find a bug, please let us know. Contact Jun Xu at 
csjunxu@comp.polyu.edu.hk or the email provided on my website at www.wangliuqing.tk.
