import os
from codefiles.visualization.visualization_algorithms import plot_line


def main() -> None:
    """
    Main function that runs the program.
    """
    folder_path = "data/output/hillclimber_data/"

    # List all files in the folder
    files = os.listdir(folder_path)

    # Filter only CSV files
    csv_files = [file for file in files if file.endswith(".csv")]

    lst_ = []
    for csv_ in csv_files:
        file_path = os.path.join(folder_path, csv_)
        with open(file_path) as file:
            lst_.append(plot_line(file))
    output_path = 'data/output/hillclimber_data/plots/'
    os.makedirs(output_path, exist_ok=True)

    index = 0
    for plots, csv_filename in zip(lst_, csv_files):
        figure_title = csv_files[index]
        for plot in plots:
            title = plot.get_axes()[0].get_title()
            figure_name = f"{figure_title}_{title}_plot.png"
            figure_path = os.path.join(
                output_path+figure_title+"/", figure_name)
            plot.savefig(figure_path)
            print(f"Saved plot: {figure_path}")
        index += 1


if __name__ == "__main__":
    main()
