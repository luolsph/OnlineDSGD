# About

Online Inference with Debiased Stochastic Gradient Descent (Online DSGD)

This repository includes source codes for online debiased stochastic gradient descent (OnlineDSGD) algorithms. These algorithms are developed for online statistical inference with high-dimensional streaming data.

Maintainers: Ruijian Han <ruijian.han@polyu.edu.hk> and Lan Luo <l.luo@rutgers.edu>

# Description: Online statistical inference of high-dimensional streaming data. 
* [datagenerator.R] generating the simulated data

* [online_LASSO_RADAR.R] function for online debiased regularization annealed epoch dual averaging (DRADAR)
* [online_Lasso_RADAR.cpp] built-in function for online DRADAR
* [offline_LASSO_RADAR.R] function for offline DRADAR (an offline counterpart of the online version)
* [offline_LASSO_RADAR.cpp] built-in function for offline DRADAR
* [online_LASSO_ASGD.R] function for online debiased stochastic gradient descent (DSGD)
* [online_Lasso_ASGD.cpp] built-in function for online DSGD
* [eval_func.R] function for offline debiased lasso and other functions for performance evaluation

* [run_main.R] execute file for simulations 

# Citation
If you use the OnlineDSGD package, please consider citing the relevant manuscript:
Ruijian Han, Lan Luo, Yuanyuan Lin and Jian Huang (2023). Online inference with debiased gradient descent. Biometrika.
