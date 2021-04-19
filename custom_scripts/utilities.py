import numpy as np

class custom_iter: # custom iterator class that allows for retrieval of current element w/out advancing
    def __init__(self, iterable):
        self.iterator = iter(iterable)
        self.current = None
    def __next__(self):
        try:
            self.current = next(self.iterator)
        except StopIteration:
            self.current = None
        finally:
            return self.current

def split_filters(string):
    UVOT_filters = ["B","U","UVW1","UVM2","UVW2","White"]
    name_idxs = custom_iter([string.index(i) for i in UVOT_filters if i in string])
    split_list = [string[name_idxs.current:next(name_idxs)] for i in range(len(UVOT_filters))]
    split_list = [item for item in split_list if len(item)>0]
    return np.unique(split_list).tolist()   