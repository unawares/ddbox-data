import pandas
from utils.patterns import Singleton


class MosesData(metaclass=Singleton):

    def __init__(self, keys_fp: str, split_fp: str) -> None:
        self.pd = pandas.merge(pandas.read_csv(keys_fp), pandas.read_csv(split_fp), how="inner", left_on="smiles", right_on="SMILES")


def moses_data():
    return MosesData('data/keys.csv', 'data/dataset_v1.csv')
