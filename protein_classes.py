class Aminoacid:

    def __init__(self, type: str) -> None:
        self._position: tuple[int] = (0, 0)
        self._type: str = type

    def add_connection(self, other: "Aminoacid") -> None:
        self._connections.add(other)

    def check_connection(self, other: "Aminoacid") -> bool:
        return other in self._connections

    def is_adjacent(self, other: "Aminoacid") -> bool:
        return abs(self.sum_position - other.sum_position) == 1

    def get_position(self) -> tuple[int]:
        return self._position

    def set_position(self, position: tuple) -> None:
        self._position = position

    def sum_position(self) -> tuple:
        return self._position[1] + self._position[2]

    def find_type(self) -> str:
        return self._type

    def type_score(self, other: "Aminoacid") -> int:
        if (self.find_type == "P" or other.find_type == "P"):
            return 0

        if (self.find_type == "H" or other.find_type == "H"):
            return -1

        if (self.find_type == "C"):
            return -5


class Grid:

    def __init__(self):
        pass
