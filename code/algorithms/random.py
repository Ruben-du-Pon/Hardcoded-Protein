import random
from typing import Optional
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid

DIRECTIONS_2D = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]

DIRECTIONS_3D = [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
              (0, -1, 0), (0, 0, 1), (0, 0, -1)]


class RandomSearch:

    def __init__(self, protein: Protein, dimensions: int, no_crossing: Optional[bool] = False) -> None:
        self._protein = protein
        self._positions = {(0, 0, 0)}
        self._no_crossing: no_crossing
        self._dimensions = dimensions

    def fold_protein(self) -> None:
        current = self._protein.get_list().link
        while current.link is not None:
            if self._no_crossing:
                while current.position in self._positions:
                    self.set_position(current)
                self._positions.add(current.position)
                current = current.link
            else:
                self.set_position(current)
                current = current.link 
                

    def get_random_direction(self) -> tuple[int, int, int]:
        if self._dimensions == 2:
          directions = DIRECTIONS_2D 
        elif self._dimensions == 3:
          directions = DIRECTIONS_3D 
        return random.choice(directions)

    def set_position(self, acid: Aminoacid):
        acid.position = tuple(
            x + y + z for x, y, z in zip(acid.predecessor.position,
                                         self.get_random_direction()))
