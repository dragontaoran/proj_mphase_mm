1. 20191126: First try. Code seems work. The results are similar among the designs.
2. 20191127: Remove the random intercept to speed up the computation. Change the wave 1 design to (10, 120, 20). Change the target precesion to 0.006.
3. 20191130: Remove the random intercept does not work. Have to include the random intercept. Main results used in the paper.
4. 20200102: MI with 10 replicates instead of 5. The spread of MI designs decreases when the number of replicates increases.  
5. 20200127: Test the adaptive sample size part of the idea() function. The code seems to be correct. The model is not very stable when we sample exclusively from the second strata.  

