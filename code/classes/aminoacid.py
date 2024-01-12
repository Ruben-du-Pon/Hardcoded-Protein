from typing import Optional


class Aminoacid:

    def __init__(self, type: str, position: tuple[int],
                 predecessor: Optional["Aminoacid"] = None,
                 link: Optional["Aminoacid"] = None) -> None:
        self.position = position
        self.predecessor = predecessor
        self.link = link
        self._type: str = type

    def get_type(self) -> str:
        return self._type

    def stability_score(self, other: "Aminoacid") -> int:
        if (self.find_type == "P" or other.find_type == "P"):
            return 0

        if (self.find_type == "H" or other.find_type == "H"):
            return -1

        if (self.find_type == "C"):
            return -5

    def __str__(self):
        return self._type
