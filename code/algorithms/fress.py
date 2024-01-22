from ..classes.protein import Protein
from .random import RandomFold
from typing import List, Tuple, Dict, Optional
import copy
import pandas as pd
import math


class FressFold:
    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
        Constructs a FressFold instance for protein folding optimization.
        
        Parameters:
        protein : Protein
            The protein instance on which the folding optimization will be performed.
        """
        self._protein = protein
        self._dimensions = dimensions
        self._verbose = verbose
        self._iterations = iterations

        self.optimization_parameters: Dict[str, float] = {}
        self.stability_scores: List[float] = []

        # keep the best protein fold
        self.best_protein_structure: Protein = protein


    def run(self) -> None:
        """
        Perform the FRESS algortim on a random protein sequence.
        """  # noqa
        best_score = 0
        best_random_score = 0
        best_protein = None
        best_random_protein = None
        for _ in range(1):
            protein_copy = self._clear_deepcopy()
            fold = RandomFold(protein_copy, self._dimensions, avoid_overlap=True, verbose=False).run()
            current_random_score = fold.get_score()
            if current_random_score < best_random_score:
                best_random_score = current_random_score
                best_random_protein = fold

        print(f"Best random score: {best_random_score}")
        print("Running _analyzeChain()")
        self._analyzeChain(best_random_protein)
        best_protein = best_random_protein # TEMP
        return best_protein

    def _analyzeChain(self, input_protein) -> None:
        """
        Analyzes the protein chain to identify regions for potential optimization.
        This version focuses on the connections of each individual amino acid.
        """
        print(f"Total score: {input_protein.get_score()}")
        chain_analysis = []
        current = input_protein._head

        while current:
            connections_info = []
            minus_one_count = 0

            connections = [pos for node in [current.predecessor, current.link]
                        if node and (pos := getattr(node, 'position', None))
                        and isinstance(pos, Tuple) and len(pos) == 3]

            adjacent_positions = [(1, 0, 0), (-1, 0, 0),
                                (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

            check_positions = [(c + adj[0], d + adj[1], e + adj[2]) for
                            (c, d, e), adj in zip([current.position] *
                                                    len(adjacent_positions),
                                                    adjacent_positions) if
                            (c + adj[0], d + adj[1], e + adj[2]) not in
                            connections]

            for pos_tuple in check_positions:
                if pos_tuple in input_protein._grid:
                    stability_score = current.get_stability_score(input_protein._grid[pos_tuple])
                    if stability_score == -1:
                        minus_one_count += 1

                    connection_detail = {
                        'adjacent_position': pos_tuple,
                        'stability_contribution': stability_score
                    }
                    connections_info.append(connection_detail)

            amino_acid_info = {
                'current_position': current.position,
                'connections': connections_info,
                'minus_one_count': minus_one_count
            }
            chain_analysis.append(amino_acid_info)

            current = current.link

        # Convert to DataFrame for easier analysis
        # Note: DataFrame will have nested lists, which may need special handling for some operations
        df_chain_analysis = pd.DataFrame(chain_analysis)

        # Identify 'good points'
        good_points = df_chain_analysis[df_chain_analysis['minus_one_count'] > 0]

        # Store the results
        self.good_points = good_points

        # If verbose, print the results
        if self._verbose:
            print("Good points in the protein sequence:")
            print(good_points)
            print(self._suggest_improvement(df_chain_analysis, input_protein, len(input_protein._sequence)))

    def _suggest_improvement(self, df_chain_analysis, input_protein, max_range_length, percentage=0.333):
        max_range_length = len(input_protein._sequence)
        range_length = math.ceil(max_range_length * percentage)

        # Calculate the average number of H's per segment
        total_h_count = input_protein._sequence.count('H')
        avg_h_per_segment = total_h_count / 3 - 0.15 * total_h_count # Since we have 3 segments: START, MIDDLE, END

        # Define the indices for each segment
        start_index = range(0, range_length)
        middle_index = range(range_length, 2 * range_length)
        end_index = range(max_range_length - range_length, max_range_length)

        # Calculate average stability contribution and count H's for each segment
        start_avg = df_chain_analysis.iloc[start_index]['minus_one_count'].mean()
        middle_avg = df_chain_analysis.iloc[middle_index]['minus_one_count'].mean()
        end_avg = df_chain_analysis.iloc[end_index]['minus_one_count'].mean()

        start_h_count = input_protein._sequence[:range_length].count('H')
        middle_h_count = input_protein._sequence[range_length:2 * range_length].count('H')
        end_h_count = input_protein._sequence[max_range_length - range_length:].count('H')

        # Analyze and suggest improvements
        min_avg = min(start_avg, middle_avg, end_avg)
        if min_avg == start_avg and start_h_count >= avg_h_per_segment:
            suggestion = "Consider making changes at the START of the protein to improve stability."
        elif min_avg == middle_avg and middle_h_count >= avg_h_per_segment:
            suggestion = "Consider making changes in the MIDDLE of the protein to improve stability."
        elif min_avg == end_avg and end_h_count >= avg_h_per_segment:
            suggestion = "Consider making changes at the END of the protein to improve stability."
        else:
            suggestion = "Improvement is less clear due to evenly distributed or low hydrophobic (H) content."

        return suggestion

    def _calculate_stability(self) -> float:
        """
        Computes the stability score of the protein's current structure based on its folding.
        
        Returns:
        float
            The calculated stability score for the current protein structure.
        """
        pass

    def _clear_deepcopy(self):
        deepcopy = copy.deepcopy(self._protein)
        current = deepcopy.get_head()
        while current is not None:
            deepcopy.remove_from_grid(current.position)
            current.position = (0, 0, 0)
            current = current.link
        return deepcopy
    
    def optimizeChain(self) -> None:
        """
        Applies optimization techniques to improve the protein's folding and stability.
        """
        pass

    def _get_best_structure(self) -> Protein:
        """
        Retrieves the most stable protein structure discovered during optimization.
        
        Returns:
        Protein
            The protein structure with the best stability score.
        """
        pass

    def _iterate_optimization(self) -> None:
        """
        Conducts a single iteration of the optimization algorithm, adjusting the protein folding as needed.
        """
        pass

    def _update_optimization_parameters(self, new_params: Dict[str, float]) -> None:
        """
        Updates the set of parameters that guide the optimization process.
        
        Parameters:
        new_params : Dict[str, float]
            A dictionary containing the new parameter values to be used in optimization.
        """
        self.optimization_parameters.update(new_params)

