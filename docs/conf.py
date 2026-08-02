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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "PyCGTOOL"
copyright = "2016, James Graham"
author = "James Graham"

# The full version, including alpha/beta/rc tags
release = "2.1.0a1"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "autoapi.extension",
    "sphinx_rtd_theme",
    "sphinx.ext.viewcode",
    "myst_parser",
]

source_suffix = [".rst", ".md"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["deflist"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
suppress_warnings = [
    # AutoAPI emits ambiguous references for the repeated module-local PathLike aliases.
    "ref.python",
    # README links such as LICENSE are valid from the repository root, not the docs tree.
    "myst.xref_missing",
]


# -- Sphinx AudoAPI options --------------------------------------------------

autoapi_type = "python"
autoapi_dirs = ["../src/pycgtool"]
