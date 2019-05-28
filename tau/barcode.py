import numpy as np
import json


class Barcode:
    def __init__(self, file):
        with open(file, 'r') as f:
            self.data = json.load(f)

    def __len__(self):
        """
        Length of the object
        """
        return len(self.data)

    def get_bd(self):
        """
        Birth, death and sizes at death for all the bars in
        the barcode diagram
        """
        _births, _deaths, _sizes = np.zeros((3, len(self)))

        for i, (key, value) in enumerate(self.data.items()):
            _births[i] = np.float(value['birth'])
            _deaths[i] = np.float(value['death'])

            foo = 0
            for ele in value.values():
                if type(ele) is list:
                    if len(ele) > foo:
                        foo += 1

            _sizes[i] = foo

        # remove the one immortal componenet
        ix = _deaths > -100

        return _births[ix], _deaths[ix], _sizes[ix]

    def get_voids_at_threshold(self, threshold, get_catalog=True):
        """
        Sizes and locations of the components that are alive
        at a given threshold
        """
        _loc = []

        # bar is a dictionary for each barcode
        for k, v in self.data.items():
            hist = v['history']

            lst = []
            for kk, vv in hist.items():
                if float(kk) > threshold:
                    lst += vv
            if len(lst) > 0:
                _loc.append(lst)

        # void sizes and centers
        vs, vc = [], []
        for lst in _loc:
            vs.append(len(lst))
            vc.append(np.mean(lst, 0))

        vs = np.array(vs)
        vc = np.array(vc)
        arg_srt = np.argsort(vs)
        vs = vs[arg_srt[::-1]]
        vc = vc[arg_srt[::-1]]

        if get_catalog:
            return vs, vc

        return _loc
