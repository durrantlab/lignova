<h1 align="center">lignova</h1>

<h4 align="center">Generate high-quality docked protein–ligand complexes at scale.</h4>

<h4 align="center" style="padding-bottom: 0.5em;"><a href="https://durrantlab.github.io/lignova">Documentation</a></h4>

<p align="center">
  <a href="https://github.com/psf/black" target="_blank">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
  </a>
  <a href="https://github.com/durrantlab/lignova/releases" target="_blank">
    <img src="https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--versions-e10079.svg" alt="License">
  </a>
  <a href="https://github.com/durrantlab/lignova/blob/main/LICENSE.md" target="_blank">
    <img src="https://img.shields.io/github/license/durrantlab/lignova" alt="License">
  </a>
</p>

LIGNOVA is an open-source, automated pipeline that pairs bioactive compounds from PubChem with high-resolution protein structures from the PDB to generate large-scale docked protein–ligand complexes using GNINA.

## Installation

Clone the [repository](https://github.com/durrantlab/lignova):

```bash
git clone git@github.com:durrantlab/lignova.git
```

Install `lignova` using `pip` after moving into the directory.

```sh
pip install .
```

This will install all dependencies and `lignova` into your current Python environment.

## Development

We use [pixi](https://pixi.sh/latest/) to manage Python environments and simplify the developer workflow.
Once you have [pixi](https://pixi.sh/latest/) installed, move into `lignova` directory (e.g., `cd lignova`) and install the  environment using the command

```bash
pixi install
```

Now you can activate the new virtual environment using

```sh
pixi shell
```

## License

Code contained in this project is released under the Apache-2.0 License as specified in [`LICENSE.md`][license].
This license grants you the freedom to use, modify, and distribute it as long as you include the original copyright notice contained in [`LICENSE.md`][license] and the following disclaimer.

> Portions of this code were incorporated and adapted with permission from [lignova](https://github.com/durrantlab/lignova) by durrantlab licensed under the [Apache-2.0 License](https://github.com/durrantlab/lignova/blob/main/LICENSE.md).

[license]: https://github.com/durrantlab/lignova/blob/main/LICENSE.md
