import matplotlib.pyplot as plt
from typing import Tuple
from ..classes.protein import Protein

def plot_3d(protein: Protein, colors: Tuple[str, str, str], filename:str) -> None:
    """
    Plot a 3D representation of the protein structure.

    The function generates a 3D scatter plot of amino acid positions in the protein structure.
    Amino acids of different types are distinguished by colors. The resulting plot is saved as an image file.

    Parameters
    ----------
    protein : Protein
        The protein structure to be visualized.
    colors : Tuple[str, str, str]
        A tuple of three colors representing different amino acid types (Hydrophobic, Polar, Cysteine).
    line_number : int
        An index used to name the output image file.
    algorithm : str
        The algorithm used for folding.

    Raises
    ------
    ValueError
        If the provided protein sequence is invalid.

    Notes
    -----
    - The function assumes that the protein structure has a double-linked list representation.
    - Amino acid types are identified by characters "H" (Hydrophobic), "P" (Polar), and "C" (Cysteine).

    Example
    -------
    colors = ("red", "blue", "green")
    protein_sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    my_protein = Protein(protein_sequence)
    plot_3d(my_protein, colors, line_number=0, algorithm="my_algorithm")
    """
    # Clear the previous plot
    plt.clf()

    colors = [color.lower() for color in colors]

    if "C" not in protein._sequence:
        curr_pos = protein.get_head()
        coordinates = [
            (curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])
        ]
        colors_ = [colors[0] if curr_pos.get_type() == "H" else colors[1]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(
                ((curr_pos.position[0],
                 curr_pos.position[1], curr_pos.position[2]))
            )
            colors_ = colors_ + [colors[0]
                                 if curr_pos.get_type() == "H" else colors[1]]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        fig = plt.figure("Protein Alignment")
        ax = fig.add_subplot(projection="3d")

        ax.scatter(x_coordinates, y_coordinates,
                   z_coordinates, s=50, marker="o", c=colors_)

        for i in range(1, len(protein)):
            ax.plot(
                [x_coordinates[i - 1], x_coordinates[i]],
                [y_coordinates[i - 1], y_coordinates[i]],
                [z_coordinates[i - 1], z_coordinates[i]],
                color="black",
            )

        x_min, y_min, z_min = min(x_coordinates), min(y_coordinates), min(z_coordinates)
        x_max, y_max, z_max = max(x_coordinates), max(y_coordinates), max(z_coordinates)

        if (x_min <= y_min) and (x_min <= z_min):
            minimum = x_min
        elif (y_min <= x_min) and (y_min <= z_min):
            minimum = y_min
        else:
            minimum = z_min
        
        if (x_max >= y_max) and (x_max >= z_max):
            maximum = x_max
        elif (y_max >= x_max) and (y_max >= z_max):
            maximum = y_max
        else:
            maximum = z_max

        ax.set_xlim((x_min - 2, x_max + 2))
        ax.set_ylim((y_min - 2, y_max + 2))
        ax.set_zlim((z_min - 2, z_max + 2))
        plt.axis("off")

        legend_labels = ["H", "P"]  # Replace with your custom characters
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color=color,
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        score_text = f"Score: {protein.get_score()}"
        ax.text(x_min - 2, y_max + 2, z_max + 2, score_text, fontsize=12.5, color='red')

        plt.savefig(filename, format='svg')
        print(f"{filename} created")

    elif "C" in protein._sequence:
        curr_pos = protein.get_head()
        coordinates = [
            (curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])
        ]
        colors_ = [
            colors[0]
            if curr_pos.get_type() == "H"
            else colors[1]
            if curr_pos.get_type() == "P"
            else colors[2]
        ]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(
                ((curr_pos.position[0],
                 curr_pos.position[1], curr_pos.position[2]))
            )
            colors_ = colors_ + [
                colors[0]
                if curr_pos.get_type() == "H"
                else colors[1]
                if curr_pos.get_type() == "P"
                else colors[2]
            ]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        fig = plt.figure("Protein Alignment")
        ax = fig.add_subplot(projection="3d")

        ax.scatter(x_coordinates, y_coordinates,
                   z_coordinates, s=50, marker="o", c=colors_)

        for i in range(1, len(protein)):
            ax.plot(
                [x_coordinates[i - 1], x_coordinates[i]],
                [y_coordinates[i - 1], y_coordinates[i]],
                [z_coordinates[i - 1], z_coordinates[i]],
                color="black",
            )

        x_min, y_min, z_min = min(x_coordinates), min(y_coordinates), min(z_coordinates)
        x_max, y_max, z_max = max(x_coordinates), max(y_coordinates), max(z_coordinates)

        if (x_min <= y_min) and (x_min <= z_min):
            minimum = x_min
        elif (y_min <= x_min) and (y_min <= z_min):
            minimum = y_min
        else:
            minimum = z_min
        
        if (x_max >= y_max) and (x_max >= z_max):
            maximum = x_max
        elif (y_max >= x_max) and (y_max >= z_max):
            maximum = y_max
        else:
            maximum = z_max

        ax.set_xlim((x_min - 2, x_max + 2))
        ax.set_ylim((y_min - 2, y_max + 2))
        ax.set_zlim((z_min - 2, z_max + 2))
        plt.axis("off")

        legend_labels = ["H", "P", "C"]  # Replace with your custom characters
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color=color,
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        score_text = f"Score: {protein.get_score()}"
        ax.text(x_min - 2, y_max + 2, z_max + 2, score_text, fontsize=12.5, color='red')


        plt.savefig(filename, format='svg')
        print(f"{filename} created")
