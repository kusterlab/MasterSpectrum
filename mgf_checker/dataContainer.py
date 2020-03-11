
from abc import ABCMeta, abstractmethod


class DataContainer(metaclass=ABCMeta):
    @abstractmethod
    def __init__(self):
        pass
