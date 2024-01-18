import os
from code.visualization.visualization_random_algo import plot_random

folder_path = "data/output/baseline/"

# List all files in the folder
files = os.listdir(folder_path)

# Filter only CSV files
csv_files = [file for file in files if file.endswith(".csv")]

lst_ = []
for csv_ in csv_files:
    file_path = os.path.join(folder_path, csv_)
    with open(file_path) as file:
        lst_.append(plot_random(file))

output_path = 'data/output/baseline_plots/'
os.makedirs(output_path, exist_ok=True)

for plots, csv_filename in zip(lst_, csv_files):
    for plot, figure_title in zip(plots, csv_files):
        title = plot.get_axes()[0].get_title()
        figure_name = f"{figure_title}_{title}_plot.png"
        figure_path = os.path.join(output_path, figure_name)
        plot.savefig(figure_path)
        print(f"Saved plot: {figure_path}")
