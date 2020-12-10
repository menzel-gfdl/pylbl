def stochastic_cloud_maps(app):
    exec(open("../tests/example-stochastic-clouds.py").read())


def setup(app):
    app.connect("builder-inited", stochastic_cloud_maps)
