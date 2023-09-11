# BFS_Load_Flow_in_active_distribution_network
Implementation of Standard Backward Forward Sweep (B/FS) method and Branch Current Based B/FS (BCBBFS) in MATLAB

# 
In order to perform load flow studies on radial distribution networks, Backward-Forward Sweep (BFS) methods are used. On the other hand, when the radial network contains DG units, some modified and improved methods should be used, like the Branch Current Based Backward Forward Sweep (BCBBFS) method.
* To see different load flow methods, read this:
https://ijsee.ctb.iau.ir/article_510070_052e005a33bbcfd2339112cf3a34071d.pdf
* The method used in my code is similar to the method mentioned at the link below; however, I've used reactive power updating for PV buses in the outer loop:
https://jesit.springeropen.com/articles/10.1186/s43067-021-00031-0 

# Problem
Implementation of the Standard Backward Forward Sweep (B/FS) method and Branch Current Based Backward Forward Sweep (BCBB/FS) for solving power flow in MATLAB, and the code should be able to get any desired balanced three-phase radial network as input by reading the network's information from an Excel file.

In the first report, the problem and solution are displayed, but the BCBBFS code needed a correction, and all results needed validation with Matpower, so they are discussed in the second report.

