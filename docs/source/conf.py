# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
project = "SpyIce"
copyright = "2024, Sneha Iyer"
author = "Sneha Iyer"

version = "1.0.0"
# The full version, including alpha/beta/rc tags
# release = "1.0.0-rc"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinxarg.ext",
    "sphinx_copybutton",
    "sphinx.ext.viewcode",
    "sphinxcontrib.pseudocode",
    # "sphinx_proof",
    "sphinx_design",
    "jupyter_sphinx",
    "myst_parser",
    "sphinx_favicon",
    "nbsphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx_automodapi.automodapi",
    "sphinx.ext.mathjax",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.graphviz",
    "sphinx.ext.intersphinx",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    # "sphinxcontrib.bibtex",
]

myst_enable_extensions = [
    "amsmath",
    "attrs_block",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "tasklist",
]

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_keyword = True
napoleon_use_rtype = True

autodoc_mock_imports = [
    "mpi4py",
    "h5py",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autosummary_generate = True
# autosummary_imported_members = True
autodoc_typehints = "description"
autodoc_class_signature = "separated"
autoclass_content = "both"
# automodapi_toctreedirnm = "code/api"
automodapi_inheritance_diagram = False
automodsumm_inherited_members = True
autodoc_member_order = "groupwise"
highlight_language = "python"
# pygments_style = "sphinx"
# pygments_dark_style = "monokai"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# The master toctree document.
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_context = {
    # ...
    "default_mode": "auto"
}
# -- Options for LaTeX output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/latex.html

latex_engine = "xelatex"
latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
    "extraclassoptions": "openany",
    "passoptionstopackages": r"""
\PassOptionsToPackage{svgnames}{xcolor}""",
    "fontpkg": r"""
""",
    "preamble": r"""
\usepackage{amsmath}
\usepackage[mathrm=sym]{unicode-math}
\usepackage[titles]{tocloft}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
""",
    "fncychap": r"\usepackage[Bjornstrup]{fncychap}",
    "printindex": r"\footnotesize\raggedright\printindex",
}
latex_toplevel_sectioning = "chapter"
latex_show_pagerefs = True
latex_theme = "manual"

latex_engine = "pdflatex"
latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
    "extraclassoptions": "openany",
    "passoptionstopackages": r"""
    \PassOptionsToPackage{svgnames}{xcolor}""",
    "preamble": r"""
    \usepackage[utf8]{inputenc}
    \usepackage[T1]{fontenc}
    \usepackage{lmodern}
    \usepackage[
    onlyrm,
    nosf,
    nott,
    nofligatures,
    ]{kpfonts}
    \usepackage[scaled=0.9]{FiraMono}
    \usepackage[titles]{tocloft}
    \cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
    \setlength{\cftchapnumwidth}{0.75cm}
    \setlength{\cftsecindent}{\cftchapnumwidth}
    \setlength{\cftsecnumwidth}{1.25cm}
""",
}
# latex_toplevel_sectioning = "chapter"
latex_show_pagerefs = True
latex_theme = "manual"
