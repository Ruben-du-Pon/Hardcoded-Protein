class Aminoacid:

    def __init__(self):
        self._position = (0, 0)

    def check_position(self, other: "Aminoacid") -> int:
        return abs(self.sum_position - other.sum_position)

    def set_position(self, position: tuple) -> None:
        self._position = position

    def sum_position(self) -> tuple:
        return self._position[1] + self._position[2]


class Protein:
    def __init__(self, string):
        self._length = len(string)


class Grid:

    def __init__(self):
        pass
