from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
from .random import RandomFold
from typing import List, Tuple, Dict, Optional
import copy
import pandas as pd
import math
import random
from operator import add, sub


class FressFold(RandomFold):
    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
        Initialize a FressFold object for protein folding optimization.

        Parameters:
        - protein (Protein): The protein to be optimized.
        - dimensions (int): The number of dimensions for the fold.
        - iterations (int): The number of iterations for the optimization process.
        - verbose (Optional[bool]): If set to True, detailed logs will be printed. Defaults to False.
        """
        self._protein = protein
        self._dimensions = dimensions
        self._verbose = verbose
        self._iterations = iterations

        self._number_of_improvements = 0

    def run(self) -> Protein:
        """
        Execute the FRESS algorithm on a random protein sequence.

        Returns:
        Protein: The best folded protein after running the algorithm.
        """
        return self._start_random_chain()

    def _start_random_chain(self) -> Protein:
        """
        Initialize a random folding chain for the protein.

        Randomly folds the protein and identifies the best starting score.
        Continues until a non-zero best random score is found.

        Returns:
        Protein: Best randomly folded protein.
        """
        best_score = 0
        best_random_score = 0
        best_protein: Optional[Protein] = None
        best_random_protein: Optional[Protein] = None
        while best_random_score == 0:
            protein_copy = self._clear_copy()
            fold = RandomFold(protein_copy, self._dimensions, avoid_overlap=True, verbose=False).run()
            current_random_score = fold.get_score()
            if current_random_score < best_random_score:
                best_random_score = current_random_score
                best_random_protein = fold

        print(f"Random start score: {best_random_score}") if self._verbose else None
        if best_random_protein is not None:
            return self._analyzeChain(best_random_protein)
        else:
            raise ValueError("No valid protein found")

    def _analyzeChain(self, input_protein: Protein) -> Protein:
        """
        Analyzes the protein chain to identify regions for potential optimization.
        This version focuses on the connections of each individual amino acid.

        Parameters:
        - input_protein (Protein): The protein to be analyzed.

        Returns:
        Protein: The protein with analysis data added for further processing.
        """
        chain_analysis = []
        current: Optional[Aminoacid] = input_protein._head

        while current:
            connections_info = []
            minus_one_count = 0

            # Collecting positions connected to the current amino acid
            connections = [pos for node in [current.predecessor, current.link]
                           if node and (pos := getattr(node, 'position', None))
                           and isinstance(pos, tuple) and len(pos) == 3]

            # Defining all adjacent positions
            adjacent_positions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

            # Filtering out adjacent positions that are not already connected
            check_positions = [(c + adj[0], d + adj[1], e + adj[2]) for
                               (c, d, e), adj in zip([current.position] * len(adjacent_positions),
                                                     adjacent_positions) if
                               (c + adj[0], d + adj[1], e + adj[2]) not in connections]

            # Evaluating each position for stability contribution
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

            # Collecting information about each amino acid
            amino_acid_info = {
                'current_position': current.position,
                'connections': connections_info,
                'minus_one_count': minus_one_count
            }
            chain_analysis.append(amino_acid_info)

            current = current.link

        # Analyzing the chain to find points of interest
        df_chain_analysis = pd.DataFrame(chain_analysis)
        good_points = df_chain_analysis[df_chain_analysis['minus_one_count'] > 0]
        self.good_points = good_points

        return self._suggest_improvement(df_chain_analysis, input_protein, len(input_protein._sequence))

    def _suggest_improvement(self, df_chain_analysis: pd.DataFrame, input_protein: Protein,
                             max_range_length: int, percentage: float = 0.333) -> Protein:
        """
        Suggests potential improvements in protein folding based on analysis.

        Parameters:
        - df_chain_analysis (pd.DataFrame): DataFrame containing chain analysis data.
        - input_protein (Protein): The protein to be improved.
        - max_range_length (int): The maximum length of the range to consider for improvement.
        - percentage (float): The percentage of total length to consider for segment analysis. Defaults to 0.333.

        Returns:
        Protein: The protein with suggested improvements.
        """
        # Calculate the maximum range length based on the protein's sequence length
        max_range_length = len(input_protein._sequence)
        
        # Calculate the actual range length based on the percentage
        range_length = math.ceil(max_range_length * percentage)

        # Calculate the average number of H's per segment
        total_h_count = input_protein._sequence.count('H')
        total_c_count = input_protein._sequence.count('C')
        total_non_p_count = total_h_count + total_c_count
        if 'C' in input_protein._sequence: # Since we have 3 segments: START, MIDDLE, END
            avg_non_p_per_segment = (total_non_p_count / 3) - 1.1
        else:
            avg_non_p_per_segment = (total_non_p_count / 3) - 0.1

        # Define the indices for each segment
        start_index = range(0, range_length)
        middle_index = range(range_length, 2 * range_length)
        end_index = range(max_range_length - range_length, max_range_length)

        # Calculate average stability contribution and count H's for each segment
        start_avg = df_chain_analysis.iloc[list(start_index)]['minus_one_count'].mean()
        middle_avg = df_chain_analysis.iloc[list(middle_index)]['minus_one_count'].mean()
        end_avg = df_chain_analysis.iloc[list(end_index)]['minus_one_count'].mean()

        start_non_p_count = input_protein._sequence[:range_length].count('H')
        middle_non_p_count = input_protein._sequence[range_length:2 * range_length].count('H')
        end_non_p_count = input_protein._sequence[max_range_length - range_length:].count('H')

        # Analyze and suggest improvements
        min_avg = min(start_avg, middle_avg, end_avg)
        
        if min_avg == start_avg and float(start_non_p_count) >= avg_non_p_per_segment:
            suggestion = "Consider making changes at the START of the protein to improve stability."
        elif min_avg == middle_avg and float(middle_non_p_count) >= avg_non_p_per_segment:
            suggestion = "Consider making changes in the MIDDLE of the protein to improve stability."
        elif min_avg == end_avg and float(end_non_p_count) >= avg_non_p_per_segment:
            suggestion = "Consider making changes at the END of the protein to improve stability."
        else:
            suggestion = " NONE "
        
        # Print the suggestion if verbose mode is enabled
        if self._verbose:
            print(suggestion)

        # Check if a valid suggestion is made and the maximum number of improvements is not reached
        if " NONE " not in suggestion and self._number_of_improvements <= 30:
            self._number_of_improvements += 1
            return self._get_refold_range(input_protein, suggestion)
        elif self._number_of_improvements == 0:
            return self._start_random_chain()

        return input_protein

    def _get_refold_range(self, input_protein: Protein, suggestion: str) -> Protein:
        """
        Determines the start and end indices for refolding based on the suggestion.

        Parameters:
        - input_protein (Protein): The protein to be refolded.
        - suggestion (str): The suggestion on where to refold.

        Returns:
        Tuple[int, int]: The start and end indices for the refolding section.
        """
        # Calculate the start and end indices for refolding based on the suggestion
        if "START" in suggestion:
            refold_range = (0, math.ceil(len(self._protein._sequence) / 3))
        elif "MIDDLE" in suggestion:
            mid = len(self._protein._sequence) / 2
            refold_range = (
                math.floor(mid - len(self._protein._sequence) / 6),
                math.ceil(mid + len(self._protein._sequence) / 6)
            )
        elif "END" in suggestion:
            refold_range = (
                math.floor(len(self._protein._sequence) * 2 / 3),
                len(self._protein._sequence) - 1
            )
        print(f"refold: {refold_range}") if self._verbose else None
        # Return the refolded protein based on the calculated range
        return self._refold_section(input_protein, refold_range, 20* len(input_protein._sequence))
        

    def _calculate_stability(self, input_protein: Protein) -> float:
        """
        Computes the stability score of the protein's current structure.

        Parameters:
        - input_protein (Protein): The protein whose stability is to be calculated.

        Returns:
        float: The calculated stability score for the current protein structure.
        """
        return input_protein.get_score()
    
    def _refold_section(self, input_protein: Protein, index_range: Tuple[int, int], num_attempts: int) -> Protein:
        """
        Randomly refolds a section of the protein and keeps the best fold.

        Parameters:
        - input_protein (Protein): The protein to be refolded.
        - index_range (Tuple[int, int]): The range of indices to be refolded.
        - num_attempts (int): The number of refolding attempts to make. Defaults to 1000.

        Returns:
        Protein: The protein structure with the best fold found in this step.
        """
        best_protein = input_protein
        best_score = self._calculate_stability(input_protein)

        start_index, end_index = index_range  # Determine the refold range based on the suggestion
        for _ in range(num_attempts):
            protein_copy = copy.deepcopy(best_protein)
            
            # Randomly refold the specified section
            fold = self._random_refold(protein_copy, start_index, end_index)

            # Check if the new fold is better
            new_score = self._calculate_stability(fold)
            if new_score <= best_score:
                print(f"The lowest score found is: {new_score}") if self._verbose else None
                best_protein = fold
                best_score = new_score

        # Return the protein with the best fold
        return self._analyzeChain(best_protein)  # Replace with self._analyzeChain(best_protein) | best_protein
                                                # For recursion | once

    def _random_refold(self, protein: Protein, start_index: int, end_index: int) -> Protein:
        """
        Performs a random refolding of a specified section of the protein.

        Parameters:
        - protein (Protein): The protein to be refolded.
        - start_index (int): The starting index of the section to refold.
        - end_index (int): The ending index of the section to refold.

        Returns:
        Protein: The protein after performing the random refold.
        """
        # Get the amino acids at the specified start and end indices
        start_acid = self._get_acid_at_index(protein, start_index)
        assert start_acid is not None, "Expected start_acid to be Aminoacid, got None"

        end_acid = self._get_acid_at_index(protein, end_index)
        assert end_acid is not None, "Expected end_acid to be Aminoacid, got None"

        # Determine the start and end positions for the refolding range
        start_position = start_acid.predecessor.position if start_acid.predecessor else (0, 0, 0)
        end_position = end_acid.link.position if end_acid.link else None

        # Extract positions for the front and back ridget parts
        index = 0
        acid: Optional[Aminoacid] = protein.get_head()
        ridget_part_front = []
        ridget_part_back = []

        while acid:
            if index <= start_index:
                ridget_part_front.append(acid.position)
            if index >= end_index:
                ridget_part_back.append(acid.position)
            index += 1
            acid = acid.link

        # Calculate directions for front and back ridget parts
        front_directions = [tuple(ridget_part_front[i + 1][j] - ridget_part_front[i][j] for j in range(3))
                            for i in range(len(ridget_part_front) - 1)]
        back_directions = [tuple(ridget_part_back[i + 1][j] - ridget_part_back[i][j] for j in range(3))
                           for i in range(len(ridget_part_back) - 1)]

        # Generate new fold for the middle section
        middle_directions = self._generate_new_fold(protein, start_index, end_index)

        # Merge directions for the entire refolded region
        merged_list = []
        merged_list.extend(front_directions)
        merged_list.extend(middle_directions)
        merged_list.extend(back_directions)

        # Calculate new positions for all atoms in the refolded region
        all_positions = [(0, 0, 0)]

        for i in range(0, len(merged_list)):
            all_positions.append(tuple(map(add, all_positions[i], merged_list[i])))

        max_iter = 0
        # Ensure uniqueness of positions to avoid overlaps
        while len(all_positions) != len(set(all_positions)):
            middle_directions = self._generate_new_fold(protein, start_index, end_index)
            merged_list = []
            merged_list.extend(front_directions)
            merged_list.extend(middle_directions)
            merged_list.extend(back_directions)
            all_positions = [(0, 0, 0)]
            

            for i in range(0, len(merged_list)):
                all_positions.append(tuple(map(add, all_positions[i], merged_list[i])))
            
            if max_iter > 10**3:# inf loop restart
                return self._start_random_chain()
            max_iter += 1

        # Create a new fold using the calculated positions
        return self._create_new_fold(protein, all_positions)

    def _generate_new_fold(self, protein: Protein, start_index: int, end_index: int) -> List[Tuple[int, ...]]:
        """
        Generates a new fold path for a specified section of the protein.

        Parameters:
        - protein (Protein): The protein to be refolded.
        - start_index (int): The starting index of the section to refold.
        - end_index (int): The ending index of the section to refold.

        Returns:
        List[Tuple[int, int, int]]: A list of new positions for the specified section.
        """
        # Create a deep copy of the protein to work with
        protein_copy = self._clear_copy()

        # Initialize acid and index for traversal
        acid: Optional[Aminoacid] = protein_copy.get_head()
        index = 0
        backtracking_bool = False

        # Move to the starting index
        while index < start_index:
            if acid is not None:
                acid = acid.link
            index += 1

        path = []
        directions = super()._get_directions()

        # Determine the number of merges required based on the start and end indices
        merges = 1 if index == 0 or end_index == len(protein_copy._sequence) - 1 else 2

        # Handle the first index if it's at the beginning of the sequence
        if index == 0 and acid is not None:
            new_position = (0, 0, 0)
            acid.position = new_position
            path.append(new_position)
            acid = acid.link
            index += 1

        # Traverse the specified section and generate a new fold path
        
        while index < (end_index + merges):
            if acid is not None and acid.predecessor:
                while directions:
                    random_direction = random.choice(directions)
                    if acid and acid.predecessor and len(acid.predecessor.position) == 3:
                        new_position = (
                            acid.predecessor.position[0] + random_direction[0],
                            acid.predecessor.position[1] + random_direction[1],
                            acid.predecessor.position[2] + random_direction[2]
                        )
                    
                        if protein_copy.is_valid_fold(new_position):
                            acid.position = new_position
                            path.append(new_position)
                            break
                        else:
                            print(random_direction)
                            print(directions)
                            # Remove direction if not valid
                            if random_direction in directions:
                                directions.remove(random_direction)
                
                if not directions:
                    print("Backtracking")
                    # Backtracking is needed
                    if len(path) > 4:
                        avoid_direction_path = (
                            path[index - 2][0] - path[index - 3][0],
                            path[index - 2][1] - path[index - 3][1],
                            path[index - 2][2] - path[index - 3][2]
                        )
                    else:
                        avoid_direction_path = None
                    if len(path) > 3:
                        avoid_direction_failed_path = (
                            path[index - 1][0] - path[index - 2][0],
                            path[index - 1][1] - path[index - 2][1],
                            path[index - 1][2] - path[index - 2][2]
                        )
                    else:
                        avoid_direction_failed_path = None
                    directions = super()._get_directions()
                    
                    if avoid_direction_path and avoid_direction_path in directions:
                        directions.remove(avoid_direction_path)
                    if avoid_direction_failed_path and avoid_direction_failed_path in directions:
                        directions.remove(avoid_direction_failed_path)
                    
                    index -= 1
                    protein_copy.remove_from_grid(acid.predecessor.position)
                    acid = acid.predecessor
                    backtracking_bool = True

            if not backtracking_bool:
                if len(random_direction) == 3:
                    x, y, z = random_direction
                    avoid_direction = (-x, -y, -z)
                    directions = super()._get_directions()

                    if avoid_direction in directions:
                        directions.remove(avoid_direction)

                    index += 1
                    if acid is not None:
                        acid = acid.link
            
            backtracking_bool = False
        
        # Calculate directions between consecutive positions in the path
        path_directions = []
        for i in range(len(path) - 1):
            path_directions.append(tuple(path[i + 1][j] - path[i][j] for j in range(3)))
        
        return path_directions

    def _create_new_fold(self, protein: Protein, path: List[Tuple[int, int, int]]) -> Protein:
        """
        Creates a new folding for the protein based on the provided path.

        Parameters:
        - protein (Protein): The protein to be refolded.
        - path (List[Tuple[int, int, int]]): The new positions for the protein's folding.

        Returns:
        Protein: The protein with the new folding applied.
        """
        # Create a deep copy of the protein to work with
        protein_copy = self._clear_copy()
        acid: Optional[Aminoacid] = protein_copy.get_head()

        index = 0 
        while acid:
            new_position = path[index]
            acid.position = new_position
            protein_copy.add_to_grid(new_position, acid)
            index += 1
            if acid is not None:
                acid = acid.link

        return protein_copy

    def _get_acid_at_index(self, protein: Protein, index: int) -> Aminoacid:
        """
        Retrieves the amino acid at a specified index in the protein.

        Parameters:
        - protein (Protein): The protein from which to retrieve the amino acid.
        - index (int): The index of the amino acid to retrieve.

        Returns:
        Aminoacid: The amino acid at the specified index.
        """
        acid: Optional[Aminoacid] = protein.get_head()
        current_index = 0
        while acid and current_index < index:
            acid = acid.link
            current_index += 1
        if acid is None:
            raise ValueError("Aminoacid not found at the specified index")
        return acid

    def _clear_copy(self) -> Protein:
        """
        Creates a deep copy of the current protein and clears its grid.

        Returns:
        Protein: A deep copy of the current protein with a cleared grid.
        """
        # Create a new protein instance
        new_protein = Protein(self._protein._sequence)  # Assuming Protein() initializes a new protein object
        
        return new_protein
