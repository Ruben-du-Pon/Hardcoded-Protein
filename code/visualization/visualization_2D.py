import matplotlib.pyplot as plt
from typing import Tuple


def plot_2d(protein, colors: Tuple[str, str, str]) -> None:
    colors = [color.lower() for color in colors]

    if "C" not in protein._sequence:
        curr_pos = protein._head
        coordinates = [(curr_pos.position[0], curr_pos.position[1])]
        colors_ = [colors[0] if curr_pos._type == "H" else colors[1]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(((curr_pos.position[0], curr_pos.position[1])))
            colors_ = colors_ + [colors[0] if curr_pos._type == "H" else colors[1]]

        x_coordinates = [x for x, _ in coordinates]
        y_coordinates = [y for _, y in coordinates]

        plt.figure("Protein Alginment")

        for x, y, color in zip(x_coordinates, y_coordinates, colors_):
            plt.scatter(x, y, color=color, marker="o")

        plt.plot(x_coordinates, y_coordinates, linestyle="-", color="black", alpha=0.1)

        plt.xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        plt.ylim((min(y_coordinates) - 2, max(y_coordinates) + 2))
        plt.axis("off")

        # Create a legend outside plt.legend
        legend_labels = ["Hydrofoob", "Polair"]
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        plt.show()

    elif "C" in protein._sequence:
        curr_pos = protein._head
        coordinates = [(curr_pos.position[0], curr_pos.position[1])]
        colors_ = [
            colors[0]
            if curr_pos._type == "H"
            else colors[1]
            if curr_pos._type == "P"
            else colors[2]
        ]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(((curr_pos.position[0], curr_pos.position[1])))
            colors_ = colors_ + [
                colors[0]
                if curr_pos._type == "H"
                else colors[1]
                if curr_pos._type == "P"
                else colors[2]
            ]

        x_coordinates = [x for x, _ in coordinates]
        y_coordinates = [y for _, y in coordinates]

        plt.figure("Protein Alginment")

        for x, y, color in zip(x_coordinates, y_coordinates, colors_):
            plt.scatter(x, y, color=color, marker="o")

        plt.plot(x_coordinates, y_coordinates, linestyle="-", color="black", alpha=0.1)

        plt.xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        plt.ylim((min(y_coordinates) - 2, max(y_coordinates) + 2))
        plt.axis("off")

        # Create a legend outside plt.legend
        legend_labels = ["Hydrofoob", "Polair", "Cysteine"]
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        plt.show()
