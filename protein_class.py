from aminoacid_class import Aminoacid


class Protein:
    def __init__(self, sequence: str) -> None:
        self._sequence = sequence
        self._head = self.create_double_linked_list()
        self._grid: dict[tuple, Aminoacid] = dict()
        self._score = 0

    @staticmethod
    def create_double_linked_list(self):
        """
        This particular method can't be called for a created object.
        It's used only once, when creating the object itself.
        After the object is created, this particular function can't
        be called on top of it.
        """
        if not self._sequence:
            return None

        head = Aminoacid(predecessor=None, link=None, type=self._sequence[0])
        current = head

        for type in self._sequence[1:]:
            new_aminoacid = Aminoacid(predecessor=current, link=None, type=type)
            current._link = new_aminoacid
            current = new_aminoacid

        return head

    @classmethod
    def get_score(self):
        return self._score


# Example usage:
protein_sequence, new_seq = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH", []
print(f"\nOriginal Protein-string: {protein_sequence}\n\n")
protein = Protein(protein_sequence)
current_node = protein._head

print("From left to right: \n")
while current_node:
    print(current_node._type, end=" ")
    if current_node._link == None:
        print("\n\n")
        break
    else:
        current_node = current_node._link

print("From right to left: \n")
while current_node:
    print(current_node._type, end=" ")
    new_seq.append(current_node._type)
    if len(protein_sequence) == len(new_seq):
        print("\n\n")
        break
    else:
        current_node = current_node._predecessor
