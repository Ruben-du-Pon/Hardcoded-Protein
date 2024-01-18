import matplotlib.pyplot as plt
import csv
import math
import numpy as np


def reading_the_csv_file(name_csv: str):
    file = open(name_csv) 
    csvreader = csv.reader(file)

    rows = []
    for row in csvreader:
        rows.append(row)

    x = 1000
    c, y, avg, amount_seq = [], [], [], math.floor(len(rows)/999)

    for _ in range(amount_seq):
        for i in rows[x-1000:x]:
            y.append(int(i[-1]))
            if i[-2] not in c:
                c.append(i[-2])
        avg.append(rows[x+1])
        x += 1004

    return c, y, avg, amount_seq

c, y, avg, amount_seq = reading_the_csv_file('3D_C.csv')

fig = plt.figure(y[0])
for count, val in enumerate(range(1000, (amount_seq+1)*1000, 1000)):
    plt.hist(y[val-1000:val], bins=np.arange(min(y[val-1000:val]), max(y[val-1000:val]), 0.5))
    plt.grid(True, alpha=0.5)  # Add grid lines
    average_value = float(avg[count][-1])
    min_val = float(min(y[val-1000:val]))
    plt.axvline(average_value, color='red', linestyle='dashed', linewidth=2, label='Average: ' + str(average_value))
    plt.axvline(min_val, color='green', linestyle='dashed', linewidth=2, label='Min: ' + str(min_val))
    plt.legend()
    plt.title(c[count])
    plt.show()