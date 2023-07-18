import os


def delete_if_exists(filepath: str):
    try:
        os.remove(filepath)
    except OSError:
        pass
