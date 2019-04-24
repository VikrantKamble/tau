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

    def get_voids_at_threshold(self, threshold):
        """
        Sizes and locations of the components that are alive
        at a given threshold
        """
        _loc = []
        _size = []

        # bar is a dictionary for each barcode
        for key, value in self.data.items():
            dist = np.inf
            pos = -1

            # remove birth and death key
            keys = [key for key in value.keys() if
                    key not in ['birth', 'death']]

            # choose the closest step from the threshold
            for i in range(0, len(keys)):
                curr_dist = np.float(keys[i]) - threshold

                if (curr_dist > 0) and (curr_dist < dist):
                    pos = i

            # if such step exists, append the mean location
            if pos != -1:
                _key = keys[pos]

                _loc.append(value[_key])
                _size.append(len(value[_key]))

        return _loc, _size
