
import matplotlib.pyplot as plt


# plot the results
def plot_bar_chart(x_arr, y_arr):
    plt.plot(x_arr, y_arr)
    plt.xlabel('INDICES')
    plt.ylabel('AVERAGE HYDROPATHY')
    plt.show()


# open file and read the sequence to a string
def read_sequence(filename):
    file = open(filename)
    lines = file.readlines()
    file.close()
    lines = lines[1:]
    for i in range(len(lines)):
        lines[i] = lines[i].strip()
    data = ''.join(lines)
    return data


def calculate_average_hydropathy(data, window_size):
    # store hydrophaty scores for each amino acid in a dictionary object
    scores = {"I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
              "A": 1.8, "G": -0.4, "T": -0.7, "W": -0.9, "S": -0.8, "Y": -1.3,
              "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5,
              "K": -3.9, "R": -4.5}
    # list for storing the average of each window
    averages = []
    # list for storing the index of the element in the middle of each window
    indices = []
    # loop for finding the average of each window
    for i in range(0, len(data), window_size):
        sum = 0
        if len(data) - i >= window_size:
            for j in range(window_size):
                sum += scores[data[i+j]]
            averages.append(sum/window_size)
            indices.append(i+int(window_size/2))
        else:
            for j in range(len(data) - i):
                sum += scores[data[i+j]]
            averages.append(sum/(len(data) - i))
            indices.append(i + int((len(data) - i)/2))
    # plot the results
    plot_bar_chart(indices, averages)


def main():
    seq = read_sequence("P08100.fasta")
    calculate_average_hydropathy(seq, 5)
    calculate_average_hydropathy(seq, 20)

zx
main()