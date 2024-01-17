import random
import copy
from collections import deque
from ..classes.protein import Protein
from .random import RandomFold


class AnnealingFold:

    def __init__(self, protein: Protein, dimensions: int) -> None:
        if dimensions not in (2, 3):
            raise ValueError(
                "Please enter dimensions as 2 or 3 for 2D or 3D folding.")

        self._protein = protein
        self._dimensions = dimensions

    def run(self) -> None:
        start_state = RandomFold(self._protein, self._dimensions, True)
        start_state.run()
        highscore = (self._protein, self._protein.get_score())
        for _ in range(1000):
            self.run_experiment()
            if self._protein.get_score() > highscore[1]:
                highscore = (self._protein, self._protein.get_score())
        self._protein = highscore[0]

    def run_experiment(self) -> None:
        temperature = -5
        protein_copy = copy.deepcopy(self._protein)
        snippet = self.get_snippet(protein_copy)
        score_old = self._protein.get_score()
        self.randomize_positions(snippet)
        score_new = protein_copy.get_score()
        chance = 2.0 ** ((score_old - score_new) / temperature)
        while chance < random.random():
            temperature = temperature * 0.75
            chance = 2.0 ** ((score_old - score_new) / temperature)
            randomise_positions(snippet)
        self._protein = protein_copy

    def get_snippet(self, protein):
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

    def randomize_positions(self, sequence):
        # Extract positions
        positions = [acid.position for acid in sequence]

        # Initialize queue with the first position and an empty path
        queue = deque([(positions[0], [])])

        # Set of visited positions
        visited = set()

        # Directions: up, down, left, right
        directions = [(0, 1, 0), (0, -1, 0), (-1, 0, 0),
                      (1, 0, 0), (0, 0, 1), (0, 0, -1)]

        while queue:
            # Dequeue a position and path
            position, path = queue.popleft()

            for direction in directions:
                # Calculate new position
                new_position = (
                    position[0] + direction[0], position[1] + direction[1], position[2] + direction[2])

                # If the new position is the same as the last position, check if all positions have been visited
                if new_position == positions[-1]:
                    new_positions = path + [new_position]
                    # If all positions have been visited, assign new positions to the sequence and return
                    if len(new_positions) == len(sequence):
                        for i, acid in enumerate(sequence):
                            acid.position = new_positions[i]
                        return
                    # Otherwise, continue with the BFS

                # If the new position is valid (within bounds and not visited), add it to the queue with the path plus the new position
                if new_position not in visited:
                    queue.append((new_position, path + [new_position]))
                    visited.add(new_position)

        # If the queue is empty and we haven't returned, there's no valid path
        raise ValueError("No valid path found")
