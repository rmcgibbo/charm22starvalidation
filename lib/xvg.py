import numpy as np

class XVGReader(object):
    def __init__(self, filename):
        self._indices = []
        self._data = []

        with open(filename) as f:
            for line in f:
                self._processLine(line.strip())

    def data(self, frame):
        return np.array(self._data[frame])

    def _processLine(self, line):
        start = line[0]
        if start == '#':
            return
        if start == '@':
            return
        
        fields = line.split()
        index = int(fields[0])
        row = map(float, fields[1:])
        self._data.append(row)
        self._indices.append(index)

if __name__ == '__main__':
    d = XVGReader('force.xvg')
    print d.data(0).reshape(-1, 3)
