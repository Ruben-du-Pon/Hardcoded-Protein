import copy
import random
from typing import List, Optional
from .random import RandomFold
from ..classes.aminoacid import Aminoacid
from ..classes.protein import Protein
from ..visualization import visualization_2D  # , visualization_3D


class HillclimberFold:
    """
    A class that folds a protein using a hillclimber algorithm.

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

    Raises
    ------
    ValueError
        If the dimensions are not 2 or 3.

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
        self._highscore = (None, None)
        self._verbose = verbose

    def run(self) -> Protein:
        """
        Runs the hillclimber algorithm.

        Returns
        -------
        Protein
            The folded protein.
        """
        # Start with a random fold
        start_state = RandomFold(self._protein, self._dimensions, True)
        protein = start_state.run()

        visualization_2D.plot_2d(
            protein, ("red", "blue", "green"), "data/output/plot/hillclimber_start.png", "png")

        protein.create_csv("data/output/csv/hillclimber_start.csv")

        # Start the highscore with the score of the random fold
        self._highscore = (protein, protein.get_score())

        # Run the algorithm for the specified number of iterations
        for _ in range(self._iterations):
            protein = self._run_experiment(protein)

        visualization_2D.plot_2d(
            protein, ("red", "blue", "green"), "data/output/plot/hillclimber_end.png", "png")

        protein.create_csv("data/output/csv/hillclimber_end.csv")

        # Return the highest scoring protein
        return self._highscore[0]

    def _run_experiment(self, protein: Protein) -> Protein:
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
        # Get a random snippet of the protein and give it an intial fold
        snippet = self._get_snippet(protein)
        changed_protein = self._fold_snippet(snippet, protein)

        # Fold the snippet until a valid fold is found
        while not self._is_valid(snippet):
            changed_protein = self._fold_snippet(snippet, changed_protein)

        # Change the protein to the new fold if it is a new highscore
        if self._check_highscore(changed_protein):
            return self._highscore[0]

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
        # Initialize snippet
        snippet = []

        # Get random length and start position
        length = random.randint(3, len(protein))
        length = min(length, 20)
        start_position = random.randint(0, len(protein) - length)
        positions = []

        current = protein.get_head()

        # Find the start position
        for _ in range(start_position):
            positions.append(current.position)
            current = current.link

        # Add the amino acids to the snippet
        for _ in range(length):
            snippet.append(current)
            positions.append(current.position)
            current = current.link

        while current:
            positions.append(current.position)
            current = current.link

        return snippet

    def _fold_snippet(self, snippet: List[Aminoacid], protein: Protein) -> Protein:
        """
        Fold a snippet until a valid fold is found.

        Parameters
        ----------
        snippet : List[Aminoacid]
            The snippet to fold.
        protein : Protein
            The protein that the snippet belongs to.

        Returns
        -------
        Protein
            The protein with the refolded snippet.
        """
        # Copy the snippet and protein
        snippet_copy = copy.deepcopy(snippet)
        protein_copy = copy.deepcopy(protein)

        # Remove the first and last amino acid from the snippet
        snippet_copy.pop(0)
        snippet_copy.pop(-1)

        # Fold the snippet
        while len(snippet_copy) > 0:
            head = snippet_copy[0]
            tail = snippet_copy[-1]

            # Store the original positions
            original_head_position = head.position
            original_tail_position = tail.position

            # Get the possible directions
            directions = [(0, 1, 0), (0, -1, 0), (-1, 0, 0),
                          (1, 0, 0), (0, 0, 1), (0, 0, -1)] if self._dimensions \
                == 3 else [(0, 1, 0), (0, -1, 0), (-1, 0, 0), (1, 0, 0)]

            # Get the possible positions for the head and tail
            head_possible_positions = [position for position in
                                       (tuple(x + y for x, y in zip(head.predecessor.position, direction))
                                        for direction in directions) if position not in
                                       protein_copy.get_grid()]

            tail_possible_positions = [position for position in
                                       (tuple(x + y for x, y in zip(tail.link.position, direction))
                                        for direction in directions) if position not in
                                       protein_copy.get_grid()]

            # Get a random position for the head and tail
            while head.position in protein_copy.get_grid() and \
                    len(head_possible_positions) > 0:
                head.position = random.choice(head_possible_positions)
                head_possible_positions.remove(head.position)

            while tail.position in protein_copy.get_grid() and \
                    len(tail_possible_positions) > 0:
                tail.position = random.choice(tail_possible_positions)
                tail_possible_positions.remove(tail.position)

            protein_copy.remove_from_grid(original_head_position)
            protein_copy.remove_from_grid(original_tail_position)

            protein_copy.add_to_grid(head.position, head)
            protein_copy.add_to_grid(tail.position, tail)

            snippet_copy.pop(0)
            if len(snippet_copy) > 0:
                snippet_copy.pop(-1)

        if protein_copy.is_valid():
            snippet = snippet_copy
            return protein_copy

        return protein

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
        if self._verbose:
            print(f"Current highscore: {self._highscore[1]}")
            print(f"Current score: {protein.get_score()}")

        if protein.get_score() < self._highscore[1]:
            self._highscore = (protein, protein.get_score())
            if self._verbose:
                print(f"New highscore found: {self._highscore[1]}")
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
