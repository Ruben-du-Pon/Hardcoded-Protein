import copy
import random
from typing import List, Optional, Tuple
from .random import RandomFold
from ..classes.aminoacid import Aminoacid
from ..classes.protein import Protein


class HillclimberFold:
    """
    A class that folds a protein using a hillclimber algorithm.

    Attributes
    ----------
    _protein : Protein
        The protein to fold.
    _dimensions : int
        The dimensions of the protein.
    _iterations : int
        The number of iterations for the hillclimber algorithm.
    _highscore : Tuple[Protein, float]
        The best protein and its score achieved during the hillclimber.
    _verbose : Optional[bool]
        Flag indicating whether to output verbose information.

    Methods
    -------
    run()
        Runs the hillclimber algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
        Initializes a HillclimberFold object.

        Parameters
        ----------
        protein : Protein
            The protein to fold.
        dimensions : int
            The dimensions of the protein.
        iterations : int
            The number of iterations for the hillclimber algorithm.
        verbose : Optional[bool], default=False
            Flag indicating whether to output verbose information.

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
        self._verbose = verbose

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

        # Start the highscore with the score of the random fold
        self._highscore = (self._protein, self._protein.get_score())

        # Copy the protein
        protein = copy.deepcopy(self._protein)

        # Run the algorithm for the specified number of iterations
        for iteration in range(self._iterations):
            if self._verbose:
                protein.create_csv(
                    f"data/output/csv/hillclimber_test_{iteration}.csv")
            protein = self._run_experiment(protein)

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

        return protein

    def _get_snippet(self, protein: Protein) -> List[Aminoacid]:
        """
        Get a random snippet of the protein.

        Parameters
        ----------
        protein : Protein
            The protein.

        Returns
        -------
        List[Aminoacid]
            A random snippet of the protein.
        """
        current = protein.get_head()
        length = len(protein)
        snippet = []

        for _ in range(length):
            snippet.append(current)
            current = current.link

        return snippet

    def _fold_snippet(self, snippet: List[Aminoacid], protein: Protein) -> None:
        """
        Fold a snippet until a valid fold is found.

        Parameters
        ----------
        snippet : List[Aminoacid]
            The snippet to fold.
        protein : Protein
            The protein.

        Returns
        -------
        None
        """
        snippet.pop(0)
        snippet.pop(-1)
        while snippet:
            head = snippet[0]
            tail = snippet[-1]

            protein.remove_from_grid(head.position)
            protein.remove_from_grid(tail.position)

            while head.position in protein.get_grid() and head is not tail:
                head.position = self._get_random_position(head, "predecessor")
            while tail.position in protein.get_grid():
                tail.position = self._get_random_position(tail, "link")

            protein.add_to_grid(head.position)
            protein.add_to_grid(tail.position)

            snippet.pop(0)
            if snippet:
                snippet.pop(-1)

    def _get_random_position(self, acid: Aminoacid, direction: str) -> Tuple[int, int, int]:
        """
        Get a random position for the amino acid.

        Parameters
        ----------
        acid : Aminoacid
            The amino acid.
        direction : str
            The direction to go (link or predecessor).

        Returns
        -------
        Tuple[int, int, int]
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
        Check if the protein is a new highscore.

        Parameters
        ----------
        protein : Protein
            The protein to check.

        Returns
        -------
        bool
            True if the protein is a new highscore, False if not.
        """
        if protein.get_score() < self._highscore[1]:
            self._highscore = (protein, protein.get_score())
            return True

        return False

    def _is_valid(self, snippet: List[Aminoacid]) -> bool:
        """
        Check if the snippet is a valid fold.

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

        for i in range(len(snippet) - 1):
            if not self._check_distance(snippet[i], snippet[i + 1]):
                return False
        return True

    def _check_distance(self, acid1: Aminoacid, acid2: Aminoacid) -> bool:
        """
        Check if the distance between two amino acids is 1.

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
