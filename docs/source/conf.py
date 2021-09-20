#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Configuration file for the Sphinx documentation builder.
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------
import os
import sys

from ambiance.__init__ import __module_name__
from ambiance.__version__ import __version__

name = __module_name__
Name = name.capitalize()

sys.path.insert(0, os.path.abspath('../../src/ambiance/'))
sys.setrecursionlimit(1500)

# -- Project information -----------------------------------------------------
project = Name
copyright = '2019, A. Dettmann'
author = 'Aaron Dettmann'

# version: The short X.Y version
# release: The full version, including alpha/beta/rc tags
# version = ''
version = __version__

# ===============
# AUTOMATE THINGS
# ===============

# Update the auto-docs
os.system('bash ./dev_doc/gen_auto_doc.sh')
os.system('python ./theory/make_model_page.py')

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    # 'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
]

# Paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# Source file parsers
# source_parsers = {
#         '.md': 'recommonmark.parser.CommonMarkParser',
#         }

# The suffix(es) of source filenames.
source_suffix = ['.rst', '.md']

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None

rst_prolog = f"""
.. |name| replace:: {Name}
.. |name_bold| replace:: **{Name}**
.. |author1| replace:: Aaron Dettmann
.. |license| replace:: *Apache-2.0*
.. _PyPI: https://pypi.org/project/ambiance/
.. _Conda: https://anaconda.org/conda-forge/ambiance
.. _pip: https://pypi.org/project/pip/
.. _NumPy: https://pypi.org/project/numpy/
.. _SciPy: https://pypi.org/project/scipy/
.. _SI units: https://en.wikipedia.org/wiki/International_System_of_Units
"""

# -- Options for HTML output -------------------------------------------------

# html_theme = 'classic'
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'canonical_url': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
}

# Paths that contain custom static files (such as style sheets) relative to this directory.
html_static_path = ['_static']

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = f'{name}doc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'a4paper',

    # The font size ('10pt', '11pt' or '12pt').
    'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    'preamble': '',

    # Latex figure (float) alignment
    'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, f'{name}.tex', f'{Name} Documentation',
     'Aaron Dettmann', 'manual'),
]

# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, f'{name}', f'{Name} Documentation',
     [author], 1)
]

# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, f'{name}', f'{Name} Documentation',
     author, f'{Name}', 'One line description of project.',
     'Miscellaneous'),
]

# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

# -- Extension configuration -------------------------------------------------
