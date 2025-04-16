This repository contains a small Python script to compute maximal cliques of
"trinomial graphs." These are graphs on the vertices 1, 2, ..., n - 1 such that
the edge `i - j` exists if and only if the resultant of `x^n - x^i + 1` and
`x^n - x^j + 1` is a signed power of 2.

## Requirements

- Python 3.8 or later
- `flint` development headers
- github version of `python-flint`

The flint development headers probably need to be installed. On Ubuntu, they
are contained in the `libflint-dev` package. On Fedora, they are in
`flint-devel`.

Fedora:

    sudo dnf install flint-devel

Ubuntu:

    sudo apt install libflint-dev

## Installation

    git clone git@github.com:rdbliss/trinomials.git
    cd trinomials

I recommend using a python virtual environment to avoid dependency conflicts.
(We need a development version of `python-flint`, which might conflict with
your system.)

    python -m venv trinomEnv
    source trinomEnv/bin/activate
    pip install -r requirements.txt

When you are using the package, execute `deactivate`. You'll need to run the
`source` command whenever you want to start using the package.

## Usage

    graph.py n

The script takes a single, positive integer argument `n`. It returns
information about the maximal cliques of the `n`th trinomial graph.

- Maximum clique: largest possible clique in the graph; that is, at least as large as any other clique.
- Maximal clique: clique which cannot be enlarged by adding a vertex
- Symmetric clique: a clique which is symmetric about `n/2`, meaning that it is
  fixed under the operation `x -> n - x` on the vertices `x`.
- Skew clique: a clique which is entirely below `n/2` or entirely above it,
  meaning `x <= n/2` for every vertex `x`.

The script determines the maximum cliques, and the maximal cliques which are
symmetric or skew. (It does *not* determine the largest possible skew or
symmetric cliques.)

Example:

    $ graph.py 10
    100%|███████████████████████████████████████████████████████████████████████████████████| 36/36 [00:00<00:00, 24260.11it/s]
    # clique number, maximal symmetric clique, maximal clique <= n/2
    5, 5, 0

    maximum cliques
    -----------
    [2, 4, 5, 6, 8]

    maximal symmetric cliques
    -----------
    [2, 4, 5, 6, 8]

    maximial skew cliques
    -----------

The output means that the largest cliques in the 10th graph have size 5. There
is only one of them, `[2, 4, 5, 6, 8]`, which happens to be symmetric. There
are no maximal cliques with vertices all <= 5.

Example:

    $ graph.py 20
    100%|█████████████████████████████████████████████████████████████████████████████████| 171/171 [00:00<00:00, 42880.90it/s]
    # clique number, largest maximal symmetric clique, largest clique <= n/2
    5, 5, 4

    maximum cliques
    -----------
    [4, 8, 10, 12, 16]
    [4, 8, 10, 11, 12]
    [8, 9, 10, 12, 16]

    maximal symmetric cliques
    -----------
    [4, 8, 10, 12, 16]

    maximial skew cliques
    -----------
    [2, 4, 5, 8]
    [4, 5, 6, 8]
    [4, 6, 7, 8]
    [4, 5, 8, 10]
    [4, 7, 8, 10]

The output means that the largest cliques in the 5th graph have size 5. There
are three cliques which achieve this, one of which is symmetric. The largest
maximal skew clique has size 4, and five such cliques attain that value.

## Outputs

The script has been precomputed for `n` up to 2500 or so. The results are saved
in tar and zip files in the `results` directory.

## Parallelism

The script runs in parallel, using as many CPU's as you have locally. This can
be quite demanding. For `n` around 1000 and bigger, you may want to run this
script on a computing cluster. (Computing `graph.py 4096` took 26 days of CPU
time!)
