import copy
import random
from typing import List, Tuple
from .random import RandomFold
from ..classes.aminoacid import Aminoacid
from ..classes.protein import Protein


class HillclimberFold():
    """
    A class that folds a protein using a hillclimber algorithm.

    Attributes
    ----------
    _protein : Protein
        The protein to fold.
    _dimensions : int
        The dimensions of the protein.

    Methods
    -------
    run()
        Runs the hillclimber algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int) -> None:
        """
        Initializes a HillclimberFold object.

        Parameters
        ----------
        protein : Protein
            The protein to fold.
        dimensions : int
            The dimensions of the protein.

        Raises
        ------
        ValueError
            If the dimensions are not 2 or 3.
        """
        if dimensions not in (2, 3):
            raise ValueError(
                "Please enter dimensions as 2 or 3 for 2D or 3D folding.")

        self._protein = protein
        self._dimensions = dimensions
        self._iterations = iterations
        self._highscore = (self._protein, self._protein.get_score())

    def run(self) -> None:
        """
        Runs the hillclimber algorithm.

        Returns
        -------
        None
        """
        # Start with a random fold
        start_state = RandomFold(self._protein, self._dimensions, True)
        start_state.run()
        # Copy the protein
        protein = copy.deepcopy(self._protein)

        # Run the algorithm for the specified number of iterations
        for _ in range(self._iterations):
            self._run_experiment(protein)

        # Set the protein to the highest scoring protein
        self._protein = self._highscore[0]

    def _run_experiment(self, protein: Protein) -> None:
        """
        Runs a single experiment of the hillclimber algorithm.

        Parameters
        ----------
        protein : Protein
            The protein to fold.

        Returns
        -------
        None
        """
        # Create a copy of the protein
        protein_copy = copy.deepcopy(protein)

        # Get a random snippet of the protein
        snippet = self._get_snippet(protein_copy)

        # Fold the snippet until a valid fold is found
        while not self._is_valid(snippet):
            self._fold_snippet(snippet, protein_copy)

        # Change the protein to the new fold if it is a new highscore
        if self._check_highscore(protein_copy):
            protein = self._highscore[0]

    def _get_snippet(self, protein: Protein) -> List[Aminoacid]:
        current = protein.get_list()
        length = random.randint(1, len(protein) // 5)
        start = random.randint(0, len(protein) - length)
        snippet = []

        for _ in range(start):
            current = current.link

        for _ in range(length):
            snippet.append(current)
            current = current.link

        return snippet

    def _fold_snippet(self, snippet: List[Aminoacid], protein: Protein) -> None:

        # Keep positions of 1st and last amino acid
        snippet.pop(0)
        snippet.pop(-1)
        while snippet:

            # Take the first and last objects from the list
            head = snippet[0]
            tail = snippet[-1]

            # Remove the original positions from the grid
            protein.remove_from_grid(head.position)
            protein.remove_from_grid(tail.position)

            # Get new positions for the amino acids
            while head.position in protein.get_grid() and head is not tail:
                head.position = self._get_random_position(head, "predecessor")
            while tail.position in protein.get_grid():
                tail.position = self._get_random_position(tail, "link")

            # Add the new positions to the grid
            protein.add_to_grid(head.position)
            protein.add_to_grid(tail.position)

            # Remove the amino acids from the list
            snippet.pop(0)
            if snippet:
                snippet.pop(-1)

    def _get_random_position(self, acid: Aminoacid, direction: str) -> Tuple[int, int, int]:
        """
        Returns a random position for the amino acid.

        Parameters
        ----------
        acid : Aminoacid
            The amino acid to get a random position for.
        direction : str
            Determines if we go forward or backward in the protein.

        Returns
        -------
        tuple
            The random position.
        """
        directions = [(0, 1, 0), (0, -1, 0), (-1, 0, 0),
                      (1, 0, 0), (0, 0, 1), (0, 0, -1)] if self._dimensions \
            == 3 else [(0, 1), (0, -1), (-1, 0), (1, 0)]
        origin = acid.link.position if direction == "link" \
            else acid.predecessor.position
        return tuple(x + y for x, y in zip(origin, random.choice(directions)))

    def _check_highscore(self, protein: Protein) -> bool:
        """
        Checks if the protein is a new highscore.

        Parameters
        ----------
        protein : Protein
            The protein to check.

        Returns
        -------
        bool
            True if the protein is a new highscore, False if not.
        """
        if protein.get_score() > self._highscore[1]:
            self._highscore = (protein, protein.get_score())
            return True
        return False

    def is_valid(self, snippet: List[Aminoacid]) -> bool:
        """
        Checks if the snippet is a valid fold.

        Parameters
        ----------
        snippet : List[Aminoacid]
            The snippet to check.

        Returns
        -------
        bool
            True if the snippet is valid, False if not.
        """
        if len(snippet) < 2:
            return False

        # Check if the snippet is a valid fold
        for i in range(len(snippet) - 1):
            if not self._check_distance(snippet[i], snippet[i + 1]):
                return False
        return True

    def _check_distance(self, acid1: Aminoacid, acid2: Aminoacid) -> bool:
        """
        Checks if the distance between two amino acids is 1.

        Parameters
        ----------
        acid1 : Aminoacid
            The first amino acid.
        acid2 : Aminoacid
            The second amino acid.

        Returns
        -------
        bool
            True if the distance is 1, False if not.
        """
        distance = 0
        for i in range(self._dimensions):
            distance += abs(acid1.position[i] - acid2.position[i])
        return distance == 1
