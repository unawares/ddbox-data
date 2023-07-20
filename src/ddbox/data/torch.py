import os
import pandas as pd
from torchvision.io import read_image
from torch.utils.data import Dataset
from ddbox.data.loaders import MosesDataLoader
from typing import List


class MosesDataset(Dataset):
    
    def __init__(self, split: str = 'train', attributes: List[str] = ['smiles'], transform=None):
        self.data_loader = MosesDataLoader(split, attributes)
        self.transform = transform

    def __len__(self):
        return self.data_loader.get_total()

    def __getitem__(self, idx):
        item = self.data_loader.get(idx, 1)
        if self.transform:
            item = self.transform(item)
        return item
