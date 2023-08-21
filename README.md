# Protein function prediction

A python package for the prediction of protein functions as GO terms to participate in the [CAFA-5](https://www.kaggle.com/competitions/cafa-5-protein-function-prediction/overview) competition. 

The development and testing of the functionalities were made using python3.8, on a linux ─ubuntu 20.4─ environment.

# Install

(optional) Create a virtual environment to install the package.

```
$ mkdir -p venv && python3 -m venv ./venv && source venv/bin/activate
```

And install it.

``` 
$ python3 -m pip install .
```

# Use

Enter a subshell in the virtual environment with dependencies available:

```
$ source venv/bin/activate
```

and run python tasks from inside the subshell.

## Usage tutorials

- [Tutorial1](analyses/01_working_with_proteins.ipynb): Working with proteins and and protein structures.
- [Tutorial2](analyses/02_get_go_terms_train_set.ipynb): Computing GO terms on the training set.
- [Tutorial3](analyses/03_get_candidate_terms_test_set.ipynb): Computing candidate GO terms for model training.
- [Tutorial4](analyses/04_get_ordered_candidate_terms_by_ia.ipynb): Order candidate terms by information accrued.
- [Tutorial5](analyses/05_create_model_for_top_scored_term.ipynb): Creating a single model and benchmarking its accuracy.
- [Tutorial6](analyses/06_run_models_massive.ipynb): Massive training of models.
- [Tutorial7](analyses/07_massive_prediction.ipynb): Massive prediction of protein functions.



