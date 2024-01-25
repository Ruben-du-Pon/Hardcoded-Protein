import copy
import random
from typing import List, Optional, Tuple
from .random import RandomFold
from .bfs import BfsFold
from ..classes.protein import Protein
from tqdm import tqdm


class HillclimberFold:
    """
    A class to perform a hillclimber fold on a protein.

    Attributes
    ----------
    _protein : Protein
        The protein to fold.
    _dimensions : int
        The number of dimensions to fold in.
    _iterations : int
        The number of iterations to run the algorithm.
    _highscore : Tuple[Protein, int]
        The highest scoring protein found.
    _verbose : bool
        Whether to print the progress of the algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
        Initializes the HillclimberFold class.

        Parameters
        ----------
        protein : Protein
            The protein to fold.
        dimensions : int
            The number of dimensions to fold in.
        iterations : int
            The number of iterations to run the algorithm.
        verbose : Optional[bool]
            Whether to print the progress of the algorithm.

        Raises
        ------
        ValueError
            If dimensions is not 2 or 3.
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
        Runs the hillclimber fold algorithm.

        Returns
        -------
        Protein
            The highest scoring protein found.
        """
        # Give feedback that the algorithm has started
        print("Starting hillclimber fold.")
        # Start with a random fold
        start_state = RandomFold(self._protein, self._dimensions, True)
        protein = start_state.run()

        # Start the highscore with the score of the random fold
        starting_fold = copy.deepcopy(protein)
        self._highscore = (starting_fold, starting_fold.get_score())

        # Run the algorithm for the specified number of iterations
        for iteration in tqdm(range(self._iterations)):
            best_fold = copy.deepcopy(self._highscore[0])
            self._run_experiment(best_fold)
            best_fold.reset_grid()

            # Write data to file
            sequence = "".join([str(acid) for acid in best_fold.get_list()])
            protein.create_csv(
                f"data/output/hillclimber_data/{sequence}_{iteration}.csv")

        # Return the highest scoring protein
        return self._highscore[0]

    def _run_experiment(self, protein: Protein) -> None:
        """
        Runs an experiment on the protein.

        Parameters
        ----------
        protein : Protein
            The protein to fold.

        Returns
        -------
        Protein
            The protein with the highest score.
        """
        # Get a random snippet of the protein
        start_position, end_position = self._get_snippet(protein)

        # Get the snippet
        snippet_acids = protein.get_list()[start_position:end_position]
        sequence = "".join([str(acid) for acid in snippet_acids])
        snippet = Protein(sequence)

        # Get the coordinates of the snippet
        start_coordinates = snippet_acids[0].position
        end_coordinates = snippet_acids[-1].position

        # Process the snippet
        args = (snippet, protein, start_coordinates,
                end_coordinates, start_position)
        score = self.process_snippet(args)

        # Update the highscore if necessary
        if score < self._highscore[1]:
            self._highscore = (protein, score)

    def _get_snippet(self, protein: Protein) -> Tuple[int, int]:
        """
        Generates a random start and end position for a snippet of the protein.

        Parameters
        ----------
        protein : Protein
            The protein to fold.

        Returns
        -------
        Tuple[int, int]
            The start and end position of the snippet.
        """
        # Get a random length
        length = random.randint(3, len(protein))

        # Set length to a maximum of 10
        length = min(length, 10)

        # Get a random start position
        start_position = random.randint(0, len(protein) - length)

        # Calculate the end position
        end_position = start_position + length

        return start_position, end_position

    def process_snippet(self, args: Tuple[Protein, Protein, Tuple[int, int],
                                          Tuple[int, int], int]) -> List[int]:
        snippet, protein, start_coordinates, end_coordinates, start_position = args

        # Create a deep copy of the protein for this process
        protein_copy = copy.deepcopy(protein)

        # Perform a breadth first search on the snippet
        search = BfsFold(protein_copy, self._dimensions)
        options = search.get_possible_foldings(
            snippet, start_coordinates, end_coordinates)

        # Try all options and return the best score
        best_score = None
        for list in options:
            for index, acid in enumerate(list):
                protein_copy.get_list()[start_position +
                                        index].position = acid.position
                protein_copy.reset_grid()

                # Check if the protein is valid
                if protein_copy.is_valid():
                    score = protein_copy.get_score()
                    if best_score is None or score < best_score:
                        best_score = score

        return best_score

    def _check_highscore(self, protein: Protein) -> bool:
        """
        Checks if the protein is a new highscore and updates the highscore if
        it is.

        Parameters
        ----------
        protein : Protein
            The protein to fold.

        Returns
        -------
        bool
            Whether the protein is a new highscore.
        """
        # Reset the grid
        protein.reset_grid()

        # Check if the protein is a new highscore
        if protein.is_valid() and protein.get_score() < self._highscore[1]:
            # Store a copy of the protein
            highscore = copy.deepcopy(protein)
            self._highscore = (highscore, highscore.get_score())
            print(
                f"New highscore found: {self._highscore[1]}") if self._verbose else None
            return True

        return False
