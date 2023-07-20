from ddbox.data.loaders import MosesDataLoader
from ddbox.api.submit import (
    submit_docking,
    submit_moses,
    submission_docking_status,
    submission_moses_status,
)


loader = MosesDataLoader()

print(len(loader.get(0, 100)))
