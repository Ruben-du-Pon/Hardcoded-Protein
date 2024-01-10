class Protein:
    def __init__(self, string):
        self._length: int = len(string)
        self._grid: dict[tuple, "Aminoacid"] = dict()

    def get_length(self) -> int:
        return self._length
