from typing import Any, Callable, Mapping


class Singleton(type):

    _instances: dict = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Registry:

    functions: Mapping[str, Callable] = {}

    def __init__(self):
        pass

    def register(self, name):
        def _func(func):
            self.functions[name] = func
            return func
        return _func

    def eval(self, name, *args, **kwargs):
        return self.functions[name](*args, **kwargs)
