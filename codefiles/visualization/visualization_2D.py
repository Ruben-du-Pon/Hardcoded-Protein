import matplotlib.pyplot as plt
from typing import Tuple
from ..classes.protein import Protein
import numpy as np


def plot_2d(protein: Protein, colors: Tuple[str, str, str], filename: str, output: str = "svg") -> None:
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

    # Convert colors to lowercase for consistency
    colors = [color.lower() for color in colors]

    if "C" not in protein._sequence:
        # If there are no Cysteine residues in the sequence
        curr_pos = protein.get_head()
        coordinates = [(curr_pos.position[0], curr_pos.position[1])]
        colors_ = [colors[0] if curr_pos.get_type() == "H" else colors[1]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(((curr_pos.position[0], curr_pos.position[1])))
            colors_ = colors_ + [colors[0]
                                 if curr_pos.get_type() == "H" else colors[1]]

        x_coordinates = [x for x, _ in coordinates]
        y_coordinates = [y for _, y in coordinates]

        # Create a new figure for the plot
        plt.figure("Protein Alignment")

        # Scatter plot for amino acid positions
        for x, y, color in zip(x_coordinates, y_coordinates, colors_):
            plt.scatter(x, y, s=50, color=color, marker="o")

        check_curr = protein.get_head()
        x_cor, y_cor = [], []

        for i in range(len(x_coordinates) - 1):

            current = check_curr.link.link
            
            while current is not None:
                diff = tuple(np.array(check_curr.position) - np.array(current.position))
                if check_curr._type == "H" and current._type == "H" and (diff == (1, 0, 0) or diff == (-1, 0, 0) or diff == (0, 1, 0) or diff == (0, -1, 0)):
                    coordinates = tuple(np.array(current.position) + np.array(tuple(value / 2 for value in diff)))
                    x_cor.append(coordinates[0])
                    y_cor.append(coordinates[1])

                current = current.link
            
            check_curr = check_curr.link

        for i, j in zip(x_cor, y_cor):
            plt.text(i, j-0.05, '*', fontsize=12, color='black', ha='center', va='center')


        # Line plot connecting amino acid positions
        plt.plot(x_coordinates, y_coordinates,
                 linestyle="-", color="black", alpha=0.7)

        # Set plot limits and turn off axis
        plt.xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        plt.ylim((min(y_coordinates) - 2, max(y_coordinates) + 2))
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

        # Add a text annotation for the score
        score_text = f"Score: {protein.get_score()}"
        plt.text(min(x_coordinates) / 2.5, max(y_coordinates) +
                 1, score_text, fontsize=12.5, color='red')
        
        plt.text(0.4, 0.3, "(0, 0)", fontsize=6, color='black', ha='right', va='top', fontdict={'fontweight': 'bold', 'style': 'italic'})

        # Save the plot as an SVG file
        if output == "png":
            plt.savefig(filename, format='png')
        else:
            plt.savefig(filename, format='svg')
        print(f"{filename} created")

    elif "C" in protein._sequence:
        # If there are Cysteine residues in the sequence
        curr_pos = protein.get_head()
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

        # Create a new figure for the plot
        plt.figure("Protein Alignment")

        # Scatter plot for amino acid positions
        for x, y, color in zip(x_coordinates, y_coordinates, colors_):
            plt.scatter(x, y, s=50, color=color, marker="o")


        check_curr = protein.get_head()
        x_cor_H, x_cor_C, y_cor_H, y_cor_C = [], [], [], []

        for i in range(len(x_coordinates) - 1):

            current = check_curr.link.link
            
            while current is not None:
                diff = tuple(np.array(check_curr.position) - np.array(current.position))
                if (((check_curr._type == "H" and current._type == "H") or (check_curr._type == "C" and current._type == "H") or (check_curr._type == "H" and current._type == "C")) and (diff == (1, 0, 0) or diff == (-1, 0, 0) or diff == (0, 1, 0) or diff == (0, -1, 0))):
                    coordinates = tuple(np.array(current.position) + np.array(tuple(value / 2 for value in diff)))
                    x_cor_H.append(coordinates[0])
                    y_cor_H.append(coordinates[1])
                elif check_curr._type == "C" and current._type == "C" and (diff == (1, 0, 0) or diff == (-1, 0, 0) or diff == (0, 1, 0) or diff == (0, -1, 0)):
                    coordinates = tuple(np.array(current.position) + np.array(tuple(value / 2 for value in diff)))
                    x_cor_C.append(coordinates[0])
                    y_cor_C.append(coordinates[1])

                current = current.link
            
            check_curr = check_curr.link

        for i, j in zip(x_cor_H, y_cor_H):
            plt.text(i, j-0.05, '*', fontsize=12, color='black', ha='center', va='center')


        for i, j in zip(x_cor_C, y_cor_C):
            plt.text(i, j-0.05, '#', fontsize=12, color='black', ha='center', va='center')



        # Line plot connecting amino acid positions
        plt.plot(x_coordinates, y_coordinates,
                 linestyle="-", color="black", alpha=0.7)

        # Set plot limits and turn off axis
        plt.xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        plt.ylim((min(y_coordinates) - 2, max(y_coordinates) + 2))
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

        # Add a text annotation for the score
        score_text = f"Score: {protein.get_score()}"
        plt.text(min(x_coordinates) / 2.5, max(y_coordinates) +
                 1, score_text, fontsize=12.5, color='red')
        
        plt.text(0.4, 0.3, "(0, 0)", fontsize=6, color='black', ha='right', va='top', fontdict={'fontweight': 'bold', 'style': 'italic'})

        # Save the plot as an SVG file
        if output == "png":
            plt.savefig(filename, format='png')
        else:
            plt.savefig(filename, format='svg')
        print(f"{filename} created")
