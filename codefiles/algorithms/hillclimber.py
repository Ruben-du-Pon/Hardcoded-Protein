import copy
import random
from typing import List, Optional, Tuple
from .random import RandomFold
from .bfs import BreadthFirstSearch
from ..classes.aminoacid import Aminoacid
from ..classes.protein import Protein
from ..visualization import visualization_2D  # , visualization_3D


class HillclimberFold:
    """
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
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
        """
        # Start with a random fold
        start_state = RandomFold(self._protein, self._dimensions, True)
        protein = start_state.run()

        visualization_2D.plot_2d(
            protein, ("red", "blue", "green"), "data/output/plot/hillclimber_start.png", "png")

        protein.create_csv("data/output/csv/hillclimber_start.csv")

        # Start the highscore with the score of the random fold
        starting_fold = copy.deepcopy(protein)
        self._highscore = (starting_fold, starting_fold.get_score())

        # Run the algorithm for the specified number of iterations
        for _ in range(self._iterations):
            print(f"Iteration {_ + 1}") if self._verbose else None
            protein = self._run_experiment(protein)

        visualization_2D.plot_2d(
            protein, ("red", "blue", "green"), "data/output/plot/hillclimber_end.png", "png")

        protein.create_csv("data/output/csv/hillclimber_end.csv")

        visualization_2D.plot_2d(
            self._highscore[0], ("red", "blue", "green"), "data/output/plot/hillclimber_highscore.png", "png")

        # Return the highest scoring protein
        return self._highscore[0]

    def _run_experiment(self, protein: Protein) -> Protein:
        """
        """
        # Get a random snippet of the protein
        start_position, end_position = self._get_snippet(protein)

        # Store a copy of the protein
        protein_copy = copy.deepcopy(protein)

        # Perform a breadth first search on the snippet
        search = BreadthFirstSearch(protein, self._dimensions)
        options = search.chunk(start_position, end_position)

        # Try all options
        for list in options:
            for index, acid in enumerate(list):
                protein.get_list()[start_position + index] = acid
                protein.reset_grid()
                self._check_highscore(protein)

        # Return the protein if it is a new highscore
        if self._check_highscore(protein):
            return protein

        # Return the original copy if a new highscore was not found
        return protein_copy

    def _get_snippet(self, protein: Protein) -> Tuple[int, int]:
        """
        """
        # Get a random length
        length = random.randint(3, len(protein))

        # Set length to a maximum of 5
        length = min(length, 5)

        # Get a random start position
        start_position = random.randint(0, len(protein) - length)

        # Calculate the end position
        end_position = start_position + length

        return start_position, end_position

    def _check_highscore(self, protein: Protein) -> bool:
        """
        """
        # Reset the grid
        protein.reset_grid()

        # Print the current highscore and score if verbose is True
        if self._verbose:
            print(f"Highscore: {self._highscore}")
            print(f"Current highscore: {self._highscore[1]}")
            print(f"Current score: {protein.get_score()}")

        # Check if the protein is a new highscore
        if protein.get_score() < self._highscore[1]:
            self._highscore = (protein, protein.get_score())
            if self._verbose:
                print(f"New highscore found: {self._highscore[1]}")
            return True

        return False
