class Aminoacid:

    def __init__(self) -> None:
        self._position: tuple[int] = (0, 0)
        self._connections: list["Aminoacid"] = {}

    def add_connection(self, other: "Aminoacid") -> None:
        self._connections.add(other)

    def check_connection(self, other: "Aminoacid") -> bool:
        return other in self._connections

    def check_position(self, other: "Aminoacid") -> int:
        return abs(self.sum_position - other.sum_position)

    def find_position(self) -> tuple[int]:
        return self._position

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
