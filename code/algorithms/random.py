import random
from typing import Optional
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid

DIRECTIONS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
              (0, -1, 0), (0, 0, 1), (0, 0, -1)]


class RandomSearch:

    def __init__(self, protein: Protein, cross_allowed: Optional[bool] = True) -> None:
        self._protein = protein
        self._positions = {(0, 0, 0)}
        self._cross_allowed: cross_allowed

    def fold_protein(self) -> None:
        current = self._protein.get_list().link
        while current.link is not None:
            if self._cross_allowed:
                self.set_position(current)
                current = current.link
            else:
                while current.position in self._positions:
                    self.set_position(current)
                self._positions.add(current.position)
                current = current.link

    def get_random_direction(self) -> tuple[int, int, int]:
        return random.choice(DIRECTIONS)

    def set_position(self, acid: Aminoacid):
        acid.position = tuple(
            x + y + z for x, y, z in zip(acid.predecessor.position,
                                         self.get_random_direction()))
