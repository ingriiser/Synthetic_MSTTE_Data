from synthcity.plugins import Plugins
from sklearn.datasets import load_diabetes
import pandas as pd
import numpy as np
import random


def gen_syn(data_path, syn_path, model, n_iter, count):
    """_summary_

    Args:
        data_path (string): path to .csv file with real data
        syn_path (string): path to the returned .csv file with synthetic dat
        model: model supported by synthcity
        n_iter (int): number of iterations in the training of the model
        count (int): number of synthetic data records

    Returns:
        .csv: csv file with count synthetic records synthesised from a model
    """
    data = pd.read_csv(data_path)

    Plugins(categories=["generic", "privacy"]).list()
    syn_model = Plugins().get(model, n_iter=n_iter)
    # it is not recommended to preprocess the data as this is done internally
    syn_model.fit(data)

    # generate synthetic data as a pandas df
    syn_data = syn_model.generate(count).dataframe()

    # replace generated patient_id
    syn_data.to_csv(syn_path, index=False)

    return syn_data


random.seed(1802)
gen_syn(
    "../data/liver_cirrhosis.csv",
    "../synth_liver_cirrhosis_tabsyncity.csv",
    "ctgan",
    100,
    20,
)
