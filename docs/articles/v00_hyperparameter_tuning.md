# 0. Hyperparameter tuning with adm

### Introduction

Using multiple algorithms and performing hyperparameters tuning is
crucial for diversifying species abundance modeling. Functions with the
*tune_abund* prefix allow users to easily perform model tuning for
different algorithms. The tuning is performed using a grid-search
approach, i.e., by evaluating model performance under many
hyperparameter combinations. To do this, the user provides a simple data
frame with hyperparameters as columns and the hyperparameter values to
be tested as rows.

In this vignette, we show how to construct and use the grid and give
details on the hyperparameters that can be tuned for each algorithm.

### Which hyperparameters can I tune?

*adm* features nine algorithms, each with its own set of
hyperparameters. The user can choose which hyperparameters to tune and
the values to test for each one. The tuning process is done with a
grid-search approach. The grid is a data frame with hyperparameters as
columns and the hyperparameter values to be tested as rows. See below
the list of algorithms and their hyperparameters.

[TABLE]

### What values the hyperparameters of my ADM can take?

#### Shallow Neural Networks

Users can tune the following hyperparameters for the Shallow Neural
Network algorithm:

- *size*: Number of units in the hidden layer. In theory, can take any
  positive integer greater or equal 1. In practice, it is limited by the
  hardware. Values too high can lead to memory issues, crashes and
  overfitting. Values too low can lead to underfitting.

- *decay*: Weight decay parameter. It can take values between 0 and 1.
  Values too high can lead to underfitting. Values too low can lead to
  overfitting.

#### Deep Neural Networks and Convolutional Neural Networks

Users can tune the following hyperparameters for the Deep and
Convolutional Neural Networks algorithm:

- *learning_rate*: The size of the gradient step taken during the
  optimization process. It can take values between 0 and 1. Values too
  high can lead to overshooting the minimum and diverging. Values too
  low can lead to slow convergence.

- *n_epochs*: Maximum number of epochs to train the model. It can take
  any positive integer greater or equal 1. Values too high can lead to
  overfitting. Values too low can lead to underfitting.

- *batch_size*: Size of the mini-batch used during the optimization
  process. It can take any positive integer greater or equal 1 and lower
  than the number of samples.

- *validation_patience*: Number of epochs with no improvement after
  which training will be stopped during validation loop. It can take any
  positive integer greater or equal 1. If the value is equal or higher
  than the number of epochs, the training will not stop.

- *fitting_patience*: Number of epochs with no improvement after which
  training will be stopped during the final model fit. It can take any
  positive integer greater or equal 1. If the value is equal or higher
  than the number of epochs, the training will not stop.

Other hyperparameters like number of hidden layers and neurons in each
layer are set with the functions *generate_arch_list*,
*generate_dnn_architecture* and *generate_cnn_architecture* (see in
“Tuning example” section below).

#### Extreme Gradient Boosting

Users can tune the following hyperparameters for the Extreme Gradient
Boosting algorithm:

- *nrounds*: Max number of boosting iterations. Can take any positive
  integer greater or equal 1. Values too high can lead to overfitting.
  Values too low can lead to underfitting.

- *max_depth*: The maximum depth of each tree. Can take any positive
  integer greater or equal 1.

- *eta*: The learning rate of the algorithm. It can take values between
  0 and 1. Values too high can lead to overshooting the minimum and
  diverging. Values too low can lead to slow convergence.

- *gamma*: Minimum loss reduction required to make a further partition
  on a leaf - node of the tree. The range depends on the loss
  (objective) chosen.

- *colsample_bytree*: Subsample ratio of columns when constructing each
  tree. It can take values between 1 and the number of predictors
  (columns)

- *min_child_weight*: Minimum sum of instance weight needed in a child.
  Can take any non-negative value.

- *subsample*: Subsample ratio of the training instance. Can take values
  between 0 and 1.

#### Generalized Additive Models

Users can tune the following hyperparameters for the Generalized
Additive Models algorithm:

- *inter*: Number of knots in x-axis. It can take any positive integer
  greater or equal 1.

#### Generalized Linear Models

Users can tune the following hyperparameters for the Generalized
Additive Models algorithm:

- *poly*: Polynomials degree for those continuous variables (i.e. used
  in predictors argument). It can take any positive integer greater or
  equal to 2

- *inter_order*: The interaction order between explanatory variables. It
  can take any positive integer greater or equal 1.

#### Generalized Boosted Regression

- *n.trees*: The total number of trees to fit. Can take any positive
  integer greater or equal 1.

- *interaction.depth*: The maximum depth of each tree. Can take any
  positive integer greater or equal 1. Values too high can lead to
  overfitting.

- *n.minobsinnode*: The minimum number of observations in the terminal
  nodes of the trees. Can take any positive integer greater or equal 1.
  Values too low or too high can lead to overfitting.

- *shrinkage*: The learning rate of the algorithm. Can take values
  between 0 and 1. Values too high can lead to overshooting the minimum
  and diverging.

#### Random Forests

Users can tune the following hyperparameters for the Random Forest
algorithm:

- *mtry*: Number of variables randomly sampled as candidates at each
  split. It can take values between 1 and the number of predictor
  variables in the dataset.

- *ntree*: Number of trees to grow. In theory, can take any positive
  integer greater or equal 1. In practice, it is limited by the
  hardware. Values too high can lead to memory issues and crashes.
  Values too low can lead to underfitting.

#### Support Vector Machines

Users can tune the following hyperparameters for the Support Vector
Machines algorithm:

- *kernel*: A string defining the kernel used in the algorithm. It can
  take any of the values described in the *kernlab* package.

- *sigma*: Either “automatic” (recommended) or the inverse kernel width
  for the Radial Basis kernel function “rbfdot” and the Laplacian kernel
  “laplacedot”.

- *C*: Cost of constraints violation. It can take any positive integer
  greater or equal 0. Values too high can lead to overfitting.

### How do I use the tuning functions with the grid? (Tuning example)

Tuning models with *adm* is straightforward. It only requires the user
to provide a grid with hyperparameters and values to be tested. To
construct the hyperparameter grid, we can use the *expand.grid* and
*list* function in a straightforward way, passing one vector of values
for each hyperparameter. Below is an example with Random Forest
hyperparameters, but the same approach can be used for all algorithms:

``` r

# Put each value of the hyperparameter in a vector, those in a list, and use expand.grid
raf_grid <- expand.grid(
  list(
    mtry = c(1, 2, 3),
    ntree = c(100, 200, 300, 400, 500)
  )
)
head(raf_grid)
#>   mtry ntree
#> 1    1   100
#> 2    2   100
#> 3    3   100
#> 4    1   200
#> 5    2   200
#> 6    3   200
```

This create a data.frame with 25 rows and 2 columns, where each row is a
combination of the hyperparameters. This object can be passed to the
*tune_abund* functions as the *grid* argument (see how below). Any grid
created this way will have the same number of rows as the product of the
number of values for each hyperparameter, e.g., if we have 3 values for
*mtry* and 4 values for *ntree*, the grid will have 3 x 4 = 12 rows.
Each row is a combination of the hyperparameters values.

To use the grid, we can pass it to the *tune_abund* function as the
*grid* argument:

``` r

library(adm)
#> Registered S3 method overwritten by 'bit':
#>   method   from  
#>   print.ri gamlss
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# Load some data
df <- adm::sppabund
df <- df %>% filter(species == "Species one")

raf_adm <- adm::tune_abund_raf(
  data = df,
  response = "ind_ha",
  metric = c("mae"),
  predictors = c("bio1", "bio12", "bio15"),
  partition = ".part1",
  grid = raf_grid # pass the grid here
)
#> Using provided grid.
#> Searching for optimal hyperparameters...
#> 
#> Fitting the best model...
#> Formula used for model fitting:
#> ind_ha ~ bio1 + bio12 + bio15
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3
#> The best model was achieved with: 
#>  mtry = 1 and ntree = 300
```

### Conclusion

Tuning in *adm* is easy. Different algorithms have distinct
hyperparameters, and the user can choose which ones to tune and the
values to test. In any case, is important to be aware of hyperparameters
limits and ranges, as well as the number of combinations that can be
generated, and their meaning. Consulting *adm* and its dependencies
documentation is a good practice to understand the algorithms and their
hyperparameters.
