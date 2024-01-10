class Protein:
    def __init__(self, string):
        self._length = len(string)

    def get_length(self) -> int:
        return self._length
