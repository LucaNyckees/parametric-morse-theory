import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (10, 10)


def plot_PD(coordinates):
    maximum = max([trio[2] for trio in coordinates])

    coord_0 = [trio for trio in coordinates if len(trio[0]) == 1]
    coord_1 = [trio for trio in coordinates if len(trio[0]) == 2]
    coord_2 = [trio for trio in coordinates if len(trio[0]) == 3]
    coord_3 = [trio for trio in coordinates if len(trio[0]) == 4]

    births_0 = [trio[1] for trio in coord_0]
    deaths_0 = [trio[2] for trio in coord_0]
    births_1 = [trio[1] for trio in coord_1]
    deaths_1 = [trio[2] for trio in coord_1]
    births_2 = [trio[1] for trio in coord_2]
    deaths_2 = [trio[2] for trio in coord_2]
    births_3 = [trio[1] for trio in coord_3]
    deaths_3 = [trio[2] for trio in coord_3]

    z = range(maximum + 1)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)

    ax1.scatter(births_0, deaths_0, color="skyblue", alpha=0.5)
    ax1.plot(z, z, color="orange")
    ax1.set_title("0-cells")
    ax2.scatter(births_1, deaths_1, color="skyblue", alpha=0.5)
    ax2.plot(z, z, color="orange")
    ax2.set_title("1-cells")
    ax3.scatter(births_2, deaths_2, color="skyblue", alpha=0.5)
    ax3.plot(z, z, color="orange")
    ax3.set_title("2-cells")
    ax4.scatter(births_3, deaths_3, color="skyblue", alpha=0.5)
    ax4.plot(z, z, color="orange")
    ax4.set_title("3-cells")
    fig.suptitle("Life coordinates of critical cells")

    for i in range(len(coord_0)):
        ax1.plot([births_0[i], births_0[i]], [births_0[i], deaths_0[i]], "m--")

    for i in range(len(coord_1)):
        ax2.plot([births_1[i], births_1[i]], [births_1[i], deaths_1[i]], "m--")

    for i in range(len(coord_2)):
        ax3.plot([births_2[i], births_2[i]], [births_2[i], deaths_2[i]], "m--")

    for i in range(len(coord_3)):
        ax4.plot([births_3[i], births_3[i]], [births_3[i], deaths_3[i]], "m--")

    plt.show()