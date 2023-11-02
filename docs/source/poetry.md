# Poetry

[Poetry](https://python-poetry.org) is a Python package manager that simplifies managing dependencies and packaging in Python.
It allows you to declare the libraries your project depends on, and it will manage (install/update) them for you.

Poetry offers several advantages over pip and conda.
For example, it provides

- a lockfile that ensures repeatable installs, meaning that identical versions of dependencies are used every time;
- a simplified configuration file that is easier to use than pip's `requirements.txt` or conda's `environment.yml` files;
- a built-in build system that makes it easy to package your project for distribution;
- an exhaustive dependency resolver that will always find a solution and provide a detailed explanation if no solution exists.

## Limitations

Unfortunately, poetry can only manage packages that are on [PyPI](https://pypi.org/).

## Installation

You can install this automatically using [`make poetry-download`](poetry-download).

## Virtual environment

We have found that setting `poetry` to keep the Python virtual environment in the project root helps unify development.
To enable this [setting](https://python-poetry.org/docs/configuration/#virtualenvsin-project), run:

```bash
poetry config virtualenvs.in-project false
```

By default, Poetry will try to use the Python version during Poetry's installation to create the virtual environment for the current project.
You can change this by following [these instructions](https://python-poetry.org/docs/managing-environments/#switching-between-environments); for example,

```bash
poetry env use python3.11
```

will create a new virtual environment using `python3.11` (as long as it is in your `PATH`).

If you want to check which virtual environment `poetry` is currently using, just run:

```bash
poetry env info
```

## Dependencies

Once you have correctly set up your virtual environment, you can use poetry to [install the same dependencies used for development](https://python-poetry.org/docs/cli/#install).

```bash
poetry install
```

This will read the `pyproject.toml` file from the current project, resolve dependencies, and install them in the virtual environment.

### add

The simplest way to add a Python package dependency is with the [`add`](https://python-poetry.org/docs/cli/#add) command.
For example, the following command would establish [NumPy](https://numpy.org/) as a dependency for the package.

```bash
poetry add numpy
```

### groups

It is common practice to separate packages into [poetry `groups`](https://python-poetry.org/docs/master/managing-dependencies/#dependency-groups).
For example, [Sphinx](https://www.sphinx-doc.org/en/master/) is commonly used to build package documentation but not the package itself.
A developer may need to have it installed, but not a downstream user.
Thus, we put this dependency into a group called `docs` as shown below.

```bash
poetry add --group docs sphinx
```

This provides additional control over Python environments and wheels.

### GitHub

Poetry can also manage bleeding-edge versions from [GitHub](https://github.com/) repositories.
Below, we show an example of installing [`MyST-NB`](https://github.com/executablebooks/MyST-NB) at a [specific commit](https://github.com/executablebooks/MyST-NB/commit/35ebd5429c1c315d2564af261b4df836651962eb).

```bash
poetry add --group docs git+https://github.com/executablebooks/MyST-NB@35ebd5429c1c315d2564af261b4df836651962eb
```

It uses the format of `git+https://github.com/<github_user>/<repo_name>@<commit_hash>`
