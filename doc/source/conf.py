from mutcleaner import __version__ as release
from pathlib import Path
import sys

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, f"{project_root}")

project = "mutcleaner"
copyright = "Xu Lab"
author = "Xu Lab"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.githubpages",
    "myst_parser",
    "sphinxcontrib.mermaid",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_autodoc_typehints",
]

# autosummary configuration
autosummary_generate = True
autosummary_generate_overwrite = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "inherited-members": False,
}

# MyST configuration
myst_enable_extensions = ["colon_fence", "deflist"]

html_theme = "pydata_sphinx_theme"
html_css_files = ["custom.css"]
html_static_path = ["_static"]
html_theme_options = {
    "show_prev_next": False,
    "use_edit_page_button": False,
    "show_nav_level": 2,
    "show_toc_level": 2,
}

myst_fence_as_directive = ["mermaid"]
