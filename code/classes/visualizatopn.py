import matplotlib.pyplot as plt
from protein import Protein

class Plot_Protein:
    def __init__(self, protein_instance: Protein):
        self.protein = protein_instance

    def plot(self, dimensionality: int, colors: tuple[str, str]):
        if dimensionality == 2:

            curr_pos = self.protein._head
            coordinates = [(curr_pos.position[0], curr_pos.position[1])]
            colors_ = ["Blue" if curr_pos._type == "H" else "Red"]

            while curr_pos.link is not None:
                curr_pos = curr_pos.link
                coordinates.append(((curr_pos.position[0], curr_pos.position[1])))
                # "H" = color[0]     "P" = color[1]
                colors_ = colors_ + [colors[0] if curr_pos._type == "H" else colors[1]]

            x_coordinates = [x for x, _ in coordinates]
            y_coordinates = [y for _, y in coordinates]

            plt.plot(x_coordinates, y_coordinates, marker="o", color=colors_)
            plt.xlim((0, self.protein._length))
            plt.ylim((0, self.protein._length))


g = Protein("HHPHHHPHPHHHPH")
plotd = Plot_Protein(g)
plotd.plot(dimensionality=2, colors=("Red", "Blue"))
