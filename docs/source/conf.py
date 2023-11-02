# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# For our custom sphinx extensions
sys.path.insert(0, os.path.abspath("./exts"))


# -- Project information -----------------------------------------------------

project = "lignova"
copyright = "2023, durrantlab"
author = "durrantlab"
html_title = "lignova"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "autoapi.extension",
    "myst_parser",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx_multiversion",
    "sphinx_design",
    "sphinxcontrib.mermaid",
    "sphinxemoji.sphinxemoji",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinx_togglebutton",
    "sphinxcontrib.bibtex",
]

suppress_warnings = ["autosectionlabel.*"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = []

# Updating master docs
root_doc = "index"

# Add mappings
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

# Include __init__ docstring for classes
autoclass_content = "both"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Including sphinx multiversion
templates_path = [
    "_templates",
]
smv_branch_whitelist = r"main"  # Only include the main branch
html_sidebars = {
    "**": [
        "sidebar/scroll-start.html",
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
        "versions.html",
    ],
}

# Manually copy over files to the root. These can then be referenced outside of the
# download directive.
# html_extra_path = []


# autoapi
autoapi_type = "python"
autoapi_generate_api_docs = True
autoapi_dirs = ["../../lignova"]
autoapi_add_toctree_entry = True
autoapi_python_class_content = "both"
autoapi_keep_files = False
autodoc_typehints = "description"

# bibtex
bibtex_bibfiles = ["refs.bib"]

# Header buttons
html_theme_options = {
    "source_repository": "https://github.com/durrantlab/lignova",
    "source_branch": "main",
    "source_directory": "docs/source/",
    "top_of_page_button": ["edit", "save", "launch"],
    "path_to_docs": "docs/source/",
    "repository_url": "https://github.com/durrantlab/lignova",
    "repository_branch": "main",
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "colab_url": "https://colab.research.google.com/",
        "deepnote_url": "https://deepnote.com/",
        "notebook_interface": "jupyterlab",
        "thebe": False,
    },
}

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "smartquotes",
    "replacements",
    "substitution",
]
