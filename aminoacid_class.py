from typing import Optional


class Aminoacid:

    def __init__(self, type: str, predecessor: Optional["Aminoacid"], link: Optional["Aminoacid"]) -> None:
        self._position: tuple[int] = (0, 0)
        self._predecessor = predecessor
        self._link = link
        self._type: str = type

    def set_position(self, position: tuple) -> None:
        self._position = position

    def get_type(self) -> str:
        return self._type

    def stability_score(self, other: "Aminoacid") -> int:
        if (self.find_type == "P" or other.find_type == "P"):
            return 0

        if (self.find_type == "H" or other.find_type == "H"):
            return -1

        if (self.find_type == "C"):
            return -5
