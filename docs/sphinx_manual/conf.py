project = '1-D Marine Particle Coagulation Model'
author = 'Burd Lab, University of Georgia'
release = 'June 2026'
extensions = ['sphinx.ext.imgmath']
html_theme = 'sphinx_rtd_theme'
html_theme_options = {'navigation_depth': 4}
html_static_path = ['_static']

# Use dvisvgm to render equations as SVG (works offline)
imgmath_image_format = 'svg'
imgmath_use_preview = True
imgmath_latex_preamble = r'\usepackage{amsmath,amssymb}'
