This folder provides CausalPath results and everything that are required to replicate the analysis. The method is available at causalpath.cs.umb.edu and it is described at

Babur, O, et al. "Causal interactions from proteomic profiles: Molecular data meet pathway knowledge." Patterns 2.6 (2021): 100257.

Results consist of multiple anaysis folders, each one containing a "parameters.txt" file (input) and multiple output files including a "results.txt" file.


To reproduce the contents of a results folder, please follow these step:

1. Go to a specific analysis folder and delete everything except the "parameters.txt" file.
2. Go to the folder that has causalpath.jar and run "java -jar causalpath.jar path/to/analysis/folder" (replace the last argument with an actual path).


To visualize result networks of a specific analysis, please follow the steps below:

1. Make sure the result folders are saved in your local storage.
2. Go to causalpath.cs.umb.edu
3. Click "View reslts from a previous analysis"
4. Select the parent analysis folder (or parent of multiple analysis folders) and upload.
5. Double click on a result network on the left tree to visualize it on the center.


If there is any issue, please consult to causalpath@cs.umb.edu.
