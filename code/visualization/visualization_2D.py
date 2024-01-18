import matplotlib.pyplot as plt
from typing import Tuple
from ..classes.protein import Protein


def plot_2d(protein: Protein, colors: Tuple[str, str, str], filename: str) -> None:
    """
    Plot a 2D representation of the protein structure.

    The function generates a 2D scatter plot of amino acid positions in the protein structure.
    Amino acids of different types are distinguished by colors. The resulting plot is saved as an image file.

    Parameters
    ----------
    protein : Protein
        The protein structure to be visualized.
    colors : Tuple[str, str, str]
        A tuple of three colors representing different amino acid types (Hydrophobic, Polar, Cysteine).
    filename : str
        The name of the output image file.

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
    plot_2d(my_protein, colors, "my_protein.png")
    """

    # Clear the previous plot
    plt.clf()

    colors = [color.lower() for color in colors]

    if "C" not in protein._sequence:
        curr_pos = protein.get_list()
        coordinates = [(curr_pos.position[0], curr_pos.position[1])]
        colors_ = [colors[0] if curr_pos.get_type() == "H" else colors[1]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(((curr_pos.position[0], curr_pos.position[1])))
            colors_ = colors_ + [colors[0]
                                 if curr_pos.get_type() == "H" else colors[1]]

        x_coordinates = [x for x, _ in coordinates]
        y_coordinates = [y for _, y in coordinates]

        plt.figure("Protein Alignment")

        for x, y, color in zip(x_coordinates, y_coordinates, colors_):
            plt.scatter(x, y, s=50, color=color, marker="o")

        plt.plot(x_coordinates, y_coordinates,
                 linestyle="-", color="black", alpha=0.7)

        x_min, y_min = min(x_coordinates), min(y_coordinates)
        x_max, y_max = max(x_coordinates), max(y_coordinates)

        if (x_min <= y_min):
            minimum = x_min
        else:
            minimum = y_min

        if (x_max >= y_max):
            maximum = x_max
        else:
            maximum = y_max

        plt.xlim((x_min - 2, x_max + 2))
        plt.ylim((y_min - 2, y_max + 2))
        plt.axis("off")

        # Create a legend outside plt.legend
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

        plt.savefig(filename)
        print(f"{filename} created")

    elif "C" in protein._sequence:
        curr_pos = protein.get_list()
        coordinates = [(curr_pos.position[0], curr_pos.position[1])]
        colors_ = [
            colors[0]
            if curr_pos.get_type() == "H"
            else colors[1]
            if curr_pos.get_type() == "P"
            else colors[2]
        ]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(((curr_pos.position[0], curr_pos.position[1])))
            colors_ = colors_ + [
                colors[0]
                if curr_pos.get_type() == "H"
                else colors[1]
                if curr_pos.get_type() == "P"
                else colors[2]
            ]

        x_coordinates = [x for x, _ in coordinates]
        y_coordinates = [y for _, y in coordinates]

        plt.figure("Protein Alignment")

        for x, y, color in zip(x_coordinates, y_coordinates, colors_):
            plt.scatter(x, y, s=50, color=color, marker="o")

        plt.plot(x_coordinates, y_coordinates,
                 linestyle="-", color="black", alpha=0.7)

        x_min, y_min = min(x_coordinates), min(y_coordinates)
        x_max, y_max = max(x_coordinates), max(y_coordinates)

        if (x_min <= y_min):
            minimum = x_min
        else:
            minimum = y_min

        if (x_max >= y_max):
            maximum = x_max
        else:
            maximum = y_max

        plt.xlim((x_min - 2, x_max + 2))
        plt.ylim((y_min - 2, y_max + 2))
        plt.axis("off")

        # Create a legend outside plt.legend
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

        plt.savefig(filename)
        print(f"{filename} created")
