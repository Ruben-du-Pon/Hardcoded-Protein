import copy
import random
from ..classes.protein import Protein
from .hillclimber import HillclimberFold


class AnnealingFold(HillclimberFold):
    """
    A class that folds a protein using a simulated annealing algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int) -> None:
        super().__init__(protein, dimensions, iterations)
        self._temperature = 1000
        self._cooling_rate = 0.003
        self._old_score = self._protein.get_score()

    def _run_experiment(self, protein: Protein) -> None:
        """
        Runs a single experiment of the simulated annealing algorithm.

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
        self._old_score = protein_copy.get_score()

        # Get a random snippet of the protein
        snippet = self._get_snippet(protein_copy)

        # Fold the snippet until a valid fold is found
        while not self._is_valid(snippet):
            self._fold_snippet(snippet, protein_copy)

        # Change the protein to the new fold if it is a new highscore
        if self._check_highscore(protein_copy):
            protein = self._highscore[0]

    def _check_highscore(self, protein: Protein) -> bool:
        old_score = self._old_score
        new_score = protein.get_score()
        temperature = self._temperature
        self._temperature *= 1 - self._cooling_rate
        chance = 2 ** ((old_score - new_score) / temperature)
        if (protein.get_score() > self._highscore[1]) and \
                (chance > random.random()):
            self._highscore = (protein, protein.get_score())
            return True
        return False
