class Aminoacid:

    def __init__(self):
        self._position = 0

    def check_position(self, other):
        return abs(self._position - other._position)

    def set_position(self, position):
        self._position = position


class Protein:
    def __init__(self, string):
        self._length = len(string)


class Grid:

    def __init__(self):
        pass
