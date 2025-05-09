from math import comb as binomial

import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import numpy as np
from graph import make_graph
from utils import is_dyadic
from flint import fmpq_poly
import networkx as nx


def draw_colors(n):
    G = make_graph(n)
    colors = nx.coloring.greedy_color(G, "independent_set")
    unique_colors = set(colors.values())
    graph_color_to_mpl_color = dict(zip(unique_colors, mpc.TABLEAU_COLORS))
    node_colors = [graph_color_to_mpl_color[colors[n]] for n in G.nodes()]

    pos = nx.shell_layout(G)
    nx.draw(
        G,
        pos,
        with_labels=True,
        node_color=node_colors,
        font_color="#333333",
    )


def make_percent_plot(n):
    def res_percent(k):
        return len(make_graph(k).edges) / binomial(k - 1, 2)

    def relprime_percent(k):
        return len(make_graph(k, cond="relprime").edges) / binomial(k - 1, 2)

    res = [res_percent(k) for k in range(3, n)]
    relprime = [relprime_percent(k) for k in range(3, n)]

    plt.figure(figsize=(2 * 4, 2 * 3))
    plt.style.use("ggplot")
    plt.plot(range(3, n), res, lw=3, label="dyadically resolving")
    plt.plot(range(3, n), relprime, lw=3, label="relatively prime")
    plt.title("Percent of dyadically resolving trinomial pairs")
    plt.xlabel("n")
    plt.legend()


def make_adj_plot(n):
    G = make_graph(n)
    m = nx.adjacency_matrix(G).toarray()
    plt.imshow(m)


def make_d_adj_plot(d):

    plt.figure()
    x = fmpq_poly([0, 1])
    g = x**d - 1

    def check_entry(n, i):
        f = x ** (n + 1) - x ** (i + 1) + 1
        return int(is_dyadic(int(abs(f.resultant(g)))))

    m = np.fromfunction(np.vectorize(check_entry), (d - 1, d - 1))
    plt.imshow(m)
