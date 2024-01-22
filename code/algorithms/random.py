from typing import List, Optional, Tuple
from operator import add, sub
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
import random


class RandomFold:
    """
    Class for performing a random folding on a protein sequence.

    Parameters
    ----------
    protein : Protein
        The protein sequence on which the random folding is performed.
    dimensions : int
        The number of dimensions for the folding (2 for 2D, 3 for 3D).
    avoid_overlap : Optional[bool], optional
        If True, ensures that no two amino acids occupy the same position,
        by repeatedly selecting a new position until a unique one is found.

    Attributes
    ----------
    _protein : Protein
        The protein sequence on which the random folding is performed.
    _dimensions : int
        The number of dimensions for the folding (2 for 2D, 3 for 3D).
    _avoid_overlap : Optional[bool]
        If True, ensures that no two amino acids occupy the same position.

    Methods
    -------
    run() -> None:
        Perform the random folding on the protein sequence.
        If avoid_overlap is True, ensures that no two amino acids occupy the
        same position by repeatedly selecting a new position until a unique one is found.

    get_random_direction() -> Tuple[int, int, int]:
        Get a random direction for the folding.
        Returns a tuple representing the random direction.

    set_position(acid: Aminoacid) -> None:
        Set the position of an amino acid based on a random direction.
        If avoid_overlap is True, ensures that no two amino acids occupy the same position.

    """  # noqa

    def __init__(self, protein: Protein, dimensions: int,
                 no_crossing: Optional[bool] = True, verbose: Optional[bool] = False) -> None:
        """
        Initialize the RandomFold object.

        Parameters
        ----------
        protein : Protein
            The protein structure to which the random folding algorithm is applied.
        dimensions : int
            Represents the dimensions of folding (2 or 3).
        avoid_overlap : Optional[bool], default: False
            If True, ensures that the folding algorithm avoids crossing paths.

        Raises
        ------
        ValueError
            If dimensions is not 2 or 3.

        Notes
        -----
        The RandomFold object is used to apply a random folding algorithm to the given protein.
        The algorithm starts from the second amino acid in the protein sequence and adds each
        amino acid to a grid in a random pattern.

        """
        if dimensions not in (2, 3):
            raise ValueError("Dimensions must be 2 or 3.")

        self._protein = protein
        self._dimensions = dimensions
        self._avoid_overlap = avoid_overlap
        self._verbose = verbose

    def run(self) -> None:
        """
        Perform the random folding on the protein sequence.

        If avoid_overlap is True, ensures that no two amino acids occupy the
        same position by repeatedly selecting a new position until a unique one is found.
        """  # noqa

        if not self._avoid_overlap:
            current = self._protein.get_head()
            self._protein.add_to_grid(current.position, current)

            move_history = []
            x = 0
            while current:
                self.set_position(current)
                current = current.link
                x += 1
        else:
            self.backtracking()

    def set_position(self, acid: Aminoacid) -> None:
        """
        Set the position of an amino acid based on a random direction.

        Parameters
        ----------
        acid : Aminoacid
            The amino acid for which the position is set.
        """  # noqa
        directions = self._get_directions()
        if acid.predecessor and acid.predecessor.position:
            direction = random.choice(directions)
            new_position = tuple(
                map(add, acid.predecessor.position, direction))
            acid.position = new_position
            self._protein.add_to_grid(new_position, acid)

    def backtracking(self, max_backtracking: int = 5000) -> None:
        """Backtracking"""
        backtrack_count = 0
        acid = self._protein.get_head()
        self._protein.add_to_grid(acid.position, acid)

        protein_path = [acid.position]
        for _ in range(len(self._protein._sequence) - 1):
            protein_path.append(None)

        directions = self._get_directions()

        print(f"{acid}[0]: {acid.position}") if self._verbose else None

        acid_index = 1
        acid = acid.link
        backtracking_bool = False

        while acid_index < len(self._protein._sequence):
            if acid.predecessor:
                while directions:
                    random_direction = random.choice(directions)
                    new_position = tuple(
                        map(add, acid.predecessor.position, random_direction))
                    if self._protein.is_valid_fold(new_position):
                        acid.position = new_position
                        self._protein.add_to_grid(new_position, acid)
                        protein_path[acid_index] = new_position
                        break
                    else:
                        # remove direction if not valid
                        if random_direction in directions:
                            directions.remove(random_direction)

                if not directions:
                    # backtracking
                    avoid_direction_path = tuple(
                        map(sub, protein_path[acid_index - 2], protein_path[acid_index - 3]))
                    avoid_direction_failed_path = tuple(
                        map(sub, protein_path[acid_index - 1], protein_path[acid_index - 2]))
                    directions = self._get_directions()
                    if avoid_direction_path in directions:
                        directions.remove(avoid_direction_path)
                    if avoid_direction_failed_path in directions:
                        directions.remove(avoid_direction_failed_path)
                    print("backtracking needed") if self._verbose else None
                    acid_index -= 1
                    backtrack_count += 1
                    if backtrack_count > max_backtracking:
                        current = self._protein.get_head()
                        self._protein.remove_from_grid(current.position)
                        current.position = (0, 0, 0)
                        current = current.link
                        while current is not None:
                            self._protein.remove_from_grid(current.position)
                            current.position = (0, 0, 0)
                            current = current.link

                        return self.backtracking(self._protein)
                    acid = acid.predecessor
                    backtracking_bool = True
                    print(
                        f"{acid}[{acid_index}]: {acid.position}") if self._verbose else None

            if not backtracking_bool:
                avoid_direction = tuple(-x for x in random_direction)
                directions = self._get_directions()
                if avoid_direction in directions:
                    directions.remove(avoid_direction)
                print(
                    f"{acid}[{acid_index}]: {acid.position}") if self._verbose else None
                acid_index += 1
                acid = acid.link
            backtracking_bool = False

        print(protein_path) if self._verbose else None

    def _get_directions(self) -> List[Tuple[int, int, int]]:
        if self._dimensions == 2:
            return [(1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0)]
        elif self._dimensions == 3:
            return [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                    (0, -1, 0), (0, 0, 1), (0, 0, -1)]
        return []

    def get_random_direction(self, directions: List[Tuple[int, int, int]]) -> Tuple[int, int, int]:
        """
        Get a random direction for the folding.

        Returns
        -------
        Tuple[int, int, int]
            A tuple representing the random direction.
        """
        return random.choice(directions)
