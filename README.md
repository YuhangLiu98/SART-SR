# SART-SR
This repository contains the official implementation of the paper: [**3-D Computed Laminography based on a Sequential Regularization**](https://ieeexplore.ieee.org/document/11018127) 
***************************************************************************
> **Abstract:** Accurate reconstruction of computed laminography (CL) is a challenging task because projections from the CL scan are incomplete, resulting in inter-slice aliasing and blurring. Hence, we establish an effective CL reconstruction model based on sequential regularization (SR) terms, which include the   norm of gradient along different directions and the truncated adaptive-weighted total variation (TAwTV). They facilitate edge-preserving diffusion and enhance the smoothness of the reconstructed object. To solve the proposed model, we introduce an alternating minimization algorithm that decomposes it into several sub-problems, which are solved by the Split-Bregman frame and gradient descent method. Compared with several iterative reconstruction methods, the experimental results demonstrate the effectiveness of the proposed method in terms of preserving edges, suppressing inter-slice aliasing, and reducing noise. 
***************************************************************************
### Algorithm
<div align=center>
<img src="https://github.com/YuhangLiu98/SART-SR/blob/main/img/SART-SR.png" width="500"/> 
</div>


-------
## Installation

1. SART-SR can be installed from source,
```shell
git clone https://github.com/YuhangLiu98/SART-SR.git
```
2. Then, [TIGRE](https://github.com/CERN/TIGRE) is required, for example,
```shell script
git clone https://github.com/CERN/TIGRE.git
```
3. Then, the three files for TAwTV need to be placed in the corresponding folders, and compilation information for TAwTV should be added to the 'compile.m' file. Refer to TIGRE's TV for guidance.

4. Lastly, place "SART-SR" code in the directory of "TIGRE":
```shell script
.\TIGRE-master\MATLAB\Demos\
```

-------

## Use
```
run `PCB_POCS_L0.m` to reconstruct the phantom.
```
-------

### RESULT  
<div align=center>
<img src="https://github.com/YuhangLiu98/SART-SR/blob/main/img/result1.png" width="800"/>   
<img src="https://github.com/YuhangLiu98/SART-SR/blob/main/img/result2.png" width="800"/>   
<img src="https://github.com/YuhangLiu98/SART-SR/blob/main/img/result3.png" width="800"/>   
<img src="https://github.com/YuhangLiu98/SART-SR/blob/main/img/result4.png" width="800"/>   
</div>

### Citation
```shell
@ARTICLE{11018127,
  author={Liu, Yuhang and Shu, Huazhong and Liu, Yi and Zhang, Pengcheng and Wang, Lei and Haigron, Pascal and Gui, Zhiguo},
  journal={IEEE Transactions on Nuclear Science}, 
  title={3-D Computed Laminography based on a Sequential Regularization}, 
  year={2025},
  doi={10.1109/TNS.2025.3574888}}
```
