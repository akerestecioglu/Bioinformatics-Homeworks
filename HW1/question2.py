
from collections import Counter
import matplotlib.pyplot as plt


def find_frequency(file_name):
    data = open(file_name + ".fasta")
    lines = data.readlines()
    lines = lines[1:]
    wf = Counter()
    for line in lines:
        for amino_acid in line:
            if amino_acid != '\n':
                wf[amino_acid] += 1
    total = sum(wf.values())
    for i in wf:
        wf[i] /= total
    out = open(file_name + "_aaCount.txt", "w")
    out.write(str(wf))
    out.close()
    data.close()
    return wf


def plot_bar_chart(wf):
    items = wf.items()
    items = sorted(items, key=lambda x: x[1], reverse = True)
    plt.bar(range(len(wf)), list([k[1] for k in items]), align='center')
    plt.xticks(range(len(wf)), list([k[0] for k in items]))
    plt.xlabel('Aminoacids')
    plt.ylabel('Frequency')
    plt.show()


amino_acid_frequency = find_frequency("P08100")
plot_bar_chart(amino_acid_frequency)
print(amino_acid_frequency)
