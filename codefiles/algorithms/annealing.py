import copy
import random
from typing import Optional
from .random import RandomFold
from ..classes.protein import Protein
from ..visualization import visualization_2D, visualization_3D
from .hillclimber import HillclimberFold


class AnnealingFold(HillclimberFold):
    """
    A class that folds a protein using a simulated annealing algorithm.

    Attributes
    ----------
    _protein : Protein
        The protein to fold.
    _dimensions : int
        The dimensions of the protein.
    _iterations : int
        The number of iterations to run the algorithm.
    _highscore : Tuple[Protein, int]
        The highest scoring protein and its score.
    _verbose : bool
        If True, print additional output.
    _temperature : float
        The initial temperature for the annealing process.
    _cooling_rate : float
        The rate at which the temperature decreases.
    _old_score : int
        The score of the protein before the current iteration.

    Methods
    -------
    run()
        Runs the simulated annealing algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
        Initializes an AnnealingFold object.

        Parameters
        ----------
        protein : Protein
            The protein to fold.
        dimensions : int
            The dimensions of the protein.
        iterations : int
            The number of iterations to run the algorithm.
        verbose : bool, optional
            If True, print additional output.

        Raises
        ------
        ValueError
            If the dimensions are not 2 or 3.
        """
        super().__init__(protein, dimensions, iterations, verbose)
        self._temperature = 10.0
        self._cooling_rate = 0.0005
        self._old_score = self._protein.get_score()

    def run(self) -> Protein:
        """
        Runs the Simulated Annealing fold algorithm.

        Returns
        -------
        Protein
            The highest scoring protein found.
        """
        # Start with a random fold
        start_state = RandomFold(self._protein, self._dimensions, True)
        protein = start_state.run()

        if self._dimensions == 2:
            visualization_2D.plot_2d(
                protein, ("red", "blue", "green"), "data/output/plot/annealing_start.png", "png")

        elif self._dimensions == 3:
            visualization_3D.plot_3d(
                protein, ("red", "blue", "green"), "data/output/plot/annealing_start.png", "png")

        protein.create_csv("data/output/csv/annealing_start.csv")

        # Start the highscore with the score of the random fold
        starting_fold = copy.deepcopy(protein)
        self._highscore = (starting_fold, starting_fold.get_score())

        # Run the algorithm for the specified number of iterations
        for _ in range(self._iterations):
            print(f"Iteration {_ + 1}") if self._verbose else None
            protein = self._run_experiment(protein)
            self._temperature *= 1 - self._cooling_rate

        # Return the highest scoring protein
        return self._highscore[0]

    def _check_highscore(self, protein: Protein) -> bool:
        """
        Checks if the current protein's score is a new highscore.

        Parameters
        ----------
        protein : Protein
            The protein to check.

        Returns
        -------
        bool
            True if the current protein's score is a new highscore, False otherwise.
        """
        # Reset the grid
        protein.reset_grid()

        # Get the old and new score
        old_score = self._old_score
        new_score = protein.get_score()

        # Get the temperature and update it
        temperature = self._temperature

        # Calculate the chance of accepting the new protein
        chance = 2 ** ((old_score - new_score) / temperature)

        # Check if the new protein is better than the old one and if it should be accepted
        if (protein.get_score() < self._highscore[1]) and \
                (chance > random.random()):
            self._highscore = (protein, protein.get_score())
            return True

        return False
