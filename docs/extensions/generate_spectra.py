def generate_spectra(app):
    exec(open("../tests/example-gas-optics.py").read())


def setup(app):
    app.connect("builder-inited", generate_spectra)
