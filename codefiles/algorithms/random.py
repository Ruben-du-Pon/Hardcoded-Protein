from typing import List, Optional, Tuple
from operator import add, sub
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
import random


class RandomFold:
    """
    Represents a random folding algorithm for proteins.

    Attributes:
    - protein (Protein): The protein to fold.
    - dimensions (int): The number of dimensions for folding (2 or 3).
    - avoid_overlap (Optional[bool]): Whether to avoid overlap in the folded protein (default is True).
    - verbose (Optional[bool]): Whether to print verbose output during the algorithm (default is False).

    Methods:
    - __init__(self, protein: Protein, dimensions: int, avoid_overlap: Optional[bool] = True, verbose: Optional[bool] = False) -> None:
        Initializes a new instance of the RandomFold class.
    - run(self) -> Protein:
        Runs the random folding algorithm and returns the folded protein.
    - set_position(self, acid: Aminoacid) -> None:
        Sets the position of an amino acid in the folded protein.
    - backtracking(self, max_backtracking: int = 5000) -> None:
        Performs backtracking in case of failed folding attempts.
    - _get_directions(self) -> List[Tuple[int, int, int]]:
        Returns a list of possible directions for folding based on the number of dimensions.
    - get_random_direction(self, directions: List[Tuple[int, int, int]]) -> Tuple[int, int, int]:
        Returns a random direction from the given list of directions.
    """  # noqa

    def __init__(self, protein: Protein, dimensions: int,
                 avoid_overlap: Optional[bool] = True,
                 verbose: Optional[bool] = False) -> None:
        """
        Initializes a new instance of the RandomFold class.

        Parameters:
        - protein (Protein): The protein to fold.
        - dimensions (int): The number of dimensions for folding (2 or 3).
        - avoid_overlap (Optional[bool]): Whether to avoid overlap in the folded protein (default is True).
        - verbose (Optional[bool]): Whether to print verbose output during the algorithm (default is False).

        Raises:
        - ValueError: If the dimensions parameter is not 2 or 3.
        """  # noqa
        if dimensions not in (2, 3):
            raise ValueError("Dimensions must be 2 or 3.")

        self._protein = protein
        self._dimensions = dimensions
        self._avoid_overlap = avoid_overlap
        self._verbose = verbose

    def run(self) -> Protein:
        """
        Runs the random folding algorithm and returns the folded protein.

        Parameters:
        - None

        Returns:
        - Protein: The folded protein.
        """
        if not self._avoid_overlap:
            current: Optional[Aminoacid] = self._protein.get_head()
            if current is not None:
                self._protein.add_to_grid(current.position, current)

            while current:
                self.set_position(current)
                if current is not None:
                    current = current.link
        else:
            self.backtracking()

        return self._protein

    def set_position(self, acid: Aminoacid) -> None:
        """
        Sets the position of an amino acid in the folded protein.

        Parameters:
        - acid (Aminoacid): The amino acid to set the position for.
        """
        directions = self._get_directions()
        if acid.predecessor and acid.predecessor.position:
            direction = random.choice(directions)
            new_position = tuple(
                map(add, acid.predecessor.position, direction))
            acid.position = new_position
            self._protein.add_to_grid(new_position, acid)

    def backtracking(self, max_backtracking: int = 5000) -> None:
        """
        Performs backtracking in case of failed folding attempts.

        Parameters:
        - max_backtracking (int): The maximum number of backtracking attempts allowed (default is 5000).
        """  # noqa
        backtrack_count = 0
        acid: Optional[Aminoacid] = self._protein.get_head()
        if acid is not None:
            self._protein.add_to_grid(acid.position, acid)
            protein_path: List[Tuple[int, int, int]] = [acid.position]

        for _ in range(len(self._protein._sequence) - 1):
            protein_path.append((0, 0, 0))

        directions = self._get_directions()

        print(
            f"{acid}[0]: {acid.position}") \
            if self._verbose and acid is not None else None

        acid_index = 1
        if acid is not None:
            acid = acid.link
        backtracking_bool = False

        while acid_index < len(self._protein._sequence) and acid is not None:
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
                        map(sub, protein_path[acid_index - 2],
                            protein_path[acid_index - 3]))
                    avoid_direction_failed_path = tuple(
                        map(sub, protein_path[acid_index - 1],
                            protein_path[acid_index - 2]))
                    directions = self._get_directions()

                    if avoid_direction_path in directions:
                        directions.remove(avoid_direction_path)

                    if avoid_direction_failed_path in directions:
                        directions.remove(avoid_direction_failed_path)

                    print("backtracking needed") if self._verbose else None
                    acid_index -= 1
                    backtrack_count += 1
                    if backtrack_count > max_backtracking:
                        current: Optional[Aminoacid] = self._protein.get_head()
                        if current is not None:
                            self._protein.remove_from_grid(current.position)
                            current.position = (0, 0, 0)
                            current = current.link
                            while current is not None:
                                self._protein.remove_from_grid(
                                    current.position)
                                current.position = (0, 0, 0)
                                current = current.link

                            return self.backtracking()
                    self._protein.remove_from_grid(acid.position)
                    acid = acid.predecessor
                    backtracking_bool = True
                    print(
                        f"{acid}[{acid_index}]: {acid.position}")\
                        if self._verbose else None

            if not backtracking_bool:
                if len(random_direction) == 3:
                    x, y, z = random_direction
                    avoid_direction = (-x, -y, -z)
                    directions = self._get_directions()
                    if avoid_direction in directions:
                        directions.remove(avoid_direction)
                    print(
                        f"{acid}[{acid_index}]: {acid.position}")\
                        if self._verbose else None
                    acid_index += 1
                    acid = acid.link
            backtracking_bool = False

        print(protein_path) if self._verbose else None

    def _get_directions(self) -> List[Tuple[int, int, int]]:
        """
        Returns a list of possible directions for folding based on the number of dimensions.

        Returns:
        - List[Tuple[int, int, int]]: A list of possible directions for folding.
        """  # noqa
        if self._dimensions == 2:
            return [(1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0)]
        elif self._dimensions == 3:
            return [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                    (0, -1, 0), (0, 0, 1), (0, 0, -1)]
        return []

    def get_random_direction(
        self, directions: List[Tuple[int, int, int]]
    ) -> Tuple[int, int, int]:
        """
        Returns a random direction from the given list of directions.

        Parameters:
        - directions (List[Tuple[int, int, int]]): The list of directions to choose from.

        Returns:
        - Tuple[int, int, int]: A random direction from the given list of directions. 
        """  # noqa
        return random.choice(directions)
