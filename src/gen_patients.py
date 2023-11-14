from synthcity.plugins import Plugins
from sklearn.datasets import load_diabetes
import pandas as pd
import numpy as np
import random


def gen_syn(data_path, return_path, model, n_iter, count, epsilon=1, times_n=1, seed=0):
    """_summary_

    Args:
        data_path (string): path to .csv file with real patients
        return_path (string): path to where the generated patients are written
        model: model supported by synthcity
        n_iter (int): number of iterations in the training of the model
        count (int): number of synthetic patients
        times_n (int): increase the size of data_path to simulate the performance
                        of a larger data set
        seed (int): random seed

    Returns:
        .csv: csv file with count synthetic patients synthesised from a model
    """
    data = pd.read_csv(data_path)

    if times_n != 1:
        data = data.loc[data.index.repeat(times_n)].reset_index(drop=True)

    Plugins(categories=["generic", "privacy"]).list()
    if model == "dpgan" or "pategan":
        # syn_model = Plugins().get(model, n_iter=n_iter, epsilon=epsilon)
        # default number of iterations
        syn_model = Plugins().get(
            model, n_iter=n_iter, epsilon=epsilon, random_state=seed
        )
    elif model == "privbayes":
        syn_model = Plugins().get(model, epsilon=epsilon)
    else:
        syn_model = Plugins().get(model, n_iter=n_iter)
    # it is not recommended to preprocess the data as this is done internally
    syn_model.fit(data)

    # generate synthetic data as a pandas df
    syn_data = syn_model.generate(count, random_state=seed).dataframe()

    # replace generated patient_id
    syn_data["patient_id"] = np.arange(1, count + 1)
    syn_data.to_csv(return_path, index=False)

    return syn_data


seed = 1202

epsilons = [0.5, 1, 5]
count = 0
for e in epsilons:
    for i in range(5):
        gen_syn(
            "../data/train_pat.csv",
            f"../data/synth_pat_dp/dp_pat_eps{e}_{i}.csv",
            "dpgan",
            500,
            391,
            epsilon=e,
            seed=seed + count,
        )
        count += 1
