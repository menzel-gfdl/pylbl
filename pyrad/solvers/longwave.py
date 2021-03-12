from numpy import exp, seterr, zeros


seterr(under="ignore", over="raise")


def planck(temperature, wavenumber, max_exp_arg=700.):
    """Calculates spectral radiance [W cm m-2] using Planck's Law.

    Args:
        temperature: Temperature [K].
        wavenumber: Wavenumber [cm-1].
        max_exp_arg: Argument passed to exp if overflow is detected.

    Returns:
        Spectral radiance [W cm m-2].
    """
    c1 = 1.1910429526245744e-8  # 2*h*c*c [W cm4 m-2].
    c2 = 1.4387773538277202  # h*c/kb [cm K].
    try:
        e = exp(c2*wavenumber/temperature)
    except FloatingPointError:
        e = exp(max_exp_arg)
    return (c1*wavenumber*wavenumber*wavenumber)/(e - 1.)


def effective_planck(center_temperature, edge_temperature, wavenumber, optical_depth):
    """Calculates an effective spectral radiance [W cm m-2] using the method described in
       doi: 10.1029/92JD01419

    Args:
        center_temperature: Temperature [K] at the center of the layer.
        edge_temperature: Temperature [K] at the edge of the layer.
        wavenumber: Wavenumber [cm-1].
        optical depth: Optical depth

    Returns:
        Spectral radiance [W cm m-2].
    """
    a = 0.193; b = 0.013  # Equation 16.
    b_center = planck(center_temperature, wavenumber)
    b_edge = planck(edge_temperature, wavenumber)
    return (b_center + (a*optical_depth + b*optical_depth*optical_depth)*b_edge) / \
           (1. + a*optical_depth + b*optical_depth*optical_depth)


def four_stream_no_scattering(center_temperature, edge_temperature, surface_temperature,
                              optical_depth, wavenumber, emissivity=1., max_exp_arg=700.):
    """Calculates fluxes for a one-dimensional atmospheric column, without scattering.

    Args:
        center_temperature: Array of temperatures [K] at layer centers.
        edge_temperature: Array of temperaures [K] at layer edges.
        surface_temperature: Surface temperature [K].
        optical_depth: Array of layer optical depths.
        wavenumber:  Wavenumber [cm-1].
        emissivitiy: Surface emissivity.
        max_exp_arg: Argument passed to exp if overflow is detected.

    Returns:
        Arrays of downward and upward fluxes [W cm m-2] at layer edges.
    """
    stream_parameters= [(-14.402613260847248, 0.07587638482015649),
                        (-3.0302159969901132, 0.676114979733751),
                        (-1.4925584280108841, 1.3726594476601073),
                        (-1.0746123148178333, 1.0169418413757783)]
    flux_down = zeros(optical_depth.size + 1)
    flux_up = zeros(optical_depth.size + 1);
    transmissivity = zeros(optical_depth.size)
    for c1, c2 in stream_parameters:
        for i in range(optical_depth.size):
            try:
                transmissivity[i] = exp(c1*optical_depth[i])
            except FloatingPointError:
                transmissivity[i] = exp(max_exp_arg)

        # Downward pass.
        beam = 0.
        for i in range(optical_depth.size):
            beam = beam*transmissivity[i] + \
                   (1. - transmissivity[i])*effective_planck(center_temperature[i],
                                                             edge_temperature[i+1],
                                                             wavenumber, optical_depth[i])
            flux_down[i+1] += c2*beam

        # Upward pass.
        beam = emissivity*planck(surface_temperature, wavenumber) + \
               (1. - emissivity)*beam
        flux_up[optical_depth.size] += c2*beam
        for i in range(optical_depth.size - 1, -1, -1):
            beam = beam*transmissivity[i] + \
                   (1. - transmissivity[i])*effective_planck(center_temperature[i],
                                                             edge_temperature[i],
                                                             wavenumber, optical_depth[i])
            flux_up[i] += c2*beam
    return flux_down, flux_up
