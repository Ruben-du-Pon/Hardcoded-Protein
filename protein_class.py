from typing import List, Optional

class Aminoacid:
    def __init__(self, predecessor: Optional["Aminoacid"], link: Optional["Aminoacid"], char: str) -> None:
        self._char = char
        self._predecessor = predecessor
        self._link = link

class Protein:
    def __init__(self, sequence: str) -> None:
        self._sequence = sequence
        self._head = self.create_double_linked_list(sequence)
        self._grid: dict[tuple, Aminoacid] = dict()
        self._score = 0  # Fix the initialization of _score

    def create_double_linked_list(self, sequence):
        if not sequence:
            return None

        head = Aminoacid(predecessor=None, link=None, char=sequence[0])
        current = head

        for char in sequence[1:]:
            new_aminoacid = Aminoacid(predecessor=current, link=None, char=char)
            current._link = new_aminoacid
            current = new_aminoacid

        return head

# Example usage:
protein_sequence = "ABCDEF"
protein = Protein(protein_sequence)

# Access the linked list from the instance
current_node = protein._head
while current_node:
    print(current_node._char, end=" ")
    current_node = current_node._link
    previous_node = current_node._predecessor
    print(previous_node._char, end=" ")
