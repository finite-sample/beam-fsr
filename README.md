## Beam Search for Forward Stepwise Regression

Classic forward stepwise regression is greedy: it adds the best variable one at a time. This works surprisingly well—but can go astray if early choices are suboptimal.

Beam search offers a simple fix. Instead of committing to a single best choice at each step, it retains the top-k candidates (beam width), exploring a broader slice of the model space.


### Summary of Results

We track test MSE at each step across datasets with varying degrees of noise, redundancy, and correlation.

On clean datasets, beam and greedy converge to the same MSE.

On noisy or highly correlated datasets, beam search often finds better models earlier.

Beam’s advantage lies in its ability to recover from early missteps—a key issue in variable selection.
