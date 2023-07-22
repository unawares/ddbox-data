import errno
import os
import shutil


def delete_if_exists(filepath: str):
    try:
        os.remove(filepath)
    except OSError:
        pass


def remove_path(path):
    if os.path.isfile(path) or os.path.islink(path):
        os.remove(path)  # remove the file
    elif os.path.isdir(path):
        shutil.rmtree(path)  # remove dir and all contains
    else:
        pass


def exists_file(path):
    return os.path.isfile(path)


def ensure_path(path):
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
