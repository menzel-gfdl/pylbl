from numpy import exp, sqrt, zeros


def adding(R_direct, R_diffuse, T_direct, T_diffuse, T_pure, direct_surface_albedo,
           diffuse_surface_albedo):
    """Calculates the shortwave upward and downward fluxes using the "adding" method.
       An overview of this method is provided in the appendix of https://doi.org/10.1029/94JD01310.*/

    Args:
        R_direct: Layer reflectivity for a direct beam (layer).
        R_diffuse: Layer reflectivity for a diffuse beam (layer).
        T_direct: Layer transmittance for a direct beam (layer).
        T_diffuse: Layer transmittance for a diffuse beam (layer).
        T_pure: Amount of incident radiation that passes through each layer without
                being absorbed or scattered (layer).
        direct_surface_albedo: Surface albedo for a direct beam.
        diffuse_surface_albedo: Surface albedo for a diffuse beam.

    Returns:
        R: Total upward reflectance at each level (level).
        T: Total downward transmittance at each level (level).
    """
    num_layers = R_direct.size
    num_levels = num_layers + 1

    R_direct_down = zeros(num_levels)  # Reflectance for a downward traveling direct beam.
    R_diffuse_down = zeros(num_levels)  # Reflectance for a downward traveling diffuse beam.
    R_diffuse_up = zeros(num_levels)  # Reflectance for an upward traveling diffuse beam.

    # For a downward traveling beam incident on a layer from above, calculate the
    # reflectance of the top of the layer.  Start at the lowest layer in the atmosphere,
    # and then build upward.
    R_direct_down[-1] = direct_surface_albedo
    R_diffuse_down[-1] = diffuse_surface_albedo
    for i in range(num_layers-1, -1, -1):
        a = T_pure[i]
        b = 1./(1. - R_diffuse[i]*R_diffuse_down[i+1])
        R_direct_down[i] = R_direct[i] + (a*R_direct_down[i+1] +
                                          (T_direct[i] - a)*R_diffuse_down[i+1])*T_diffuse[i]*b
        R_diffuse_down[i] = R_diffuse[i] + T_diffuse[i]*T_diffuse[i]*R_diffuse_down[i+1]*b

    # For an upward traveling beam incident on a layer from below, calculate the
    # reflectance of the bottom of the layer.  Start at the top layer in the atmosphere,
    # and then build downward.
    R_diffuse_up[0] = R_diffuse[0]
    for i in range(1, num_layers):
        a = 1./(1. - R_diffuse[i]*R_diffuse_up[i-1])
        R_diffuse_up[i] = R_diffuse[i] + T_diffuse[i]*T_diffuse[i]*R_diffuse_up[i-1]*a

    # Calculate the flux through each pressure level.
    direct_beam = 1.
    R, T = zeros(num_levels), zeros(num_levels)
    R[0] = direct_beam*R_direct_down[0]
    T[0] = direct_beam
    diffuse_beam = direct_beam*(T_direct[0] - T_pure[0])
    for i in range(1, num_levels):
        if i > 1:
            a = 1./(1. - R_diffuse[i-1]*R_diffuse_up[i-2])
            diffuse_beam = (direct_beam*R_direct[i-1]*R_diffuse_up[i-2] + diffuse_beam) * \
                T_diffuse[i-1]*a + direct_beam*(T_direct[i-1] - T_pure[i-1])
        direct_beam *= T_pure[i-1]
        b = 1./(1. - R_diffuse_down[i]*R_diffuse_up[i-1])
        R[i] = (direct_beam*R_direct_down[i] + diffuse_beam*R_diffuse_down[i])*b
        T[i] = direct_beam*(1. + R_direct_down[i]*R_diffuse_up[i-1]*b) + diffuse_beam*b
    return R, T


def delta_eddington(optical_depth, single_scatter_albedo, asymmetry_factor, cosine_zenith):
    """Uses the delta-Eddington approximation to calculate the reflectivity and transmittance
       of a plane-parallel atmospheric layer.  The specificform of the Eddinton
       approximation used here is described in
       https://doi.org/10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2.

    Args:
        optical_depth: Optical depth.
        single_scatter_albedo: Single-scatter albedo.
        asymmetry_factor: Asymmetry factor.
        cosine_zenith: Cosine of the beam zenith angle.

    Returns:
        R: Reflectivity.
        T: Transmittance.
        T_pure: Amount of incident radiation that passes through without being
                absorbed or scattered.
    """
    tau, omega, g = scaling(optical_depth, single_scatter_albedo, asymmetry_factor)
    return eddington(tau, omega, g, cosine_zenith)


def eddington(optical_depth, single_scatter_albedo, asymmetry_factor, cosine_zenith):
    """Uses the Eddington approximation to calculate the reflectivity and transmittance
       of a plane-parallel atmospheric layer.  The specificform of the Eddinton
       approximation used here is described in
       https://doi.org/10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2.

    Args:
        optical_depth: Optical depth.
        single_scatter_albedo: Single-scatter albedo.
        asymmetry_factor: Asymmetry factor.
        cosine_zenith: Cosine of the beam zenith angle.

    Returns:
        R: Reflectivity.
        T: Transmittance.
        T_pure: Amount of incident radiation that passes through without being
                absorbed or scattered.
    """
    if optical_depth <= 0.:
        # There is no gas in the layer, so no scattering or absorption.
        return 0., 1., 1.

    if single_scatter_albedo <= 0.:
        # There is no scattering, only absorption.
        R = 0.
        T = exp(-1.*optical_depth/cosine_zenith)
        T_pure = T
        return R, T, T_pure

    # Calculate the Eddington approximation parameters.
    gamma1 = 0.25*(7. - single_scatter_albedo*(4. + 3.*asymmetry_factor))  # Table 1, row 1.
    gamma2 = -0.25*(1. - single_scatter_albedo*(4. - 3.*asymmetry_factor))  # Table 1, row 1.
    gamma3 = 0.25*(2. - 3.*asymmetry_factor*cosine_zenith)  # Table 1, row 1.
    gamma4 = 1. - gamma3  # Equation 21.
    T_pure = exp(-1.*optical_depth/cosine_zenith)

    if single_scatter_albedo >= 1.:
        # Scattering is conservative, so there is no absorption.
        R = (1./(1. + gamma1*optical_depth))*(gamma1*optical_depth + (gamma3 -
                                              gamma1*cosine_zenith)*(1. - T_pure))  # Equation 24.
        T = 1. - R  # Equation 24.
        return R, T, T_pure

    alpha1 = gamma1*gamma4 + gamma2*gamma3  # Equation 16.
    alpha2 = gamma1*gamma3 + gamma2*gamma4  # Equation 17.
    k = sqrt(gamma1*gamma1 - gamma2*gamma2)  # Equation 18.

    # Prevent overflow of exponentials.
    max_exp_arg = 80.
    tau = optical_depth
    if 1./cosine_zenith > k and optical_depth/cosine_zenith > max_exp_arg:
        tau = max_exp_arg*cosine_zenith
    elif optical_depth*k > max_exp_arg:
        tau = max_exp_arg/k
    tp = exp(tau/cosine_zenith)
    tm = exp(-1.*tau/cosine_zenith)
    tkp = exp(tau*k)
    tkm = exp(-1.*tau*k)

    R = (single_scatter_albedo/((1. - k*k*cosine_zenith*cosine_zenith)*((k + gamma1)*tkp +
                                (k - gamma1)*tkm))) * \
        ((1. - k*cosine_zenith)*(alpha2 + k*gamma3)*tkp -
         (1. + k*cosine_zenith)*(alpha2 - k*gamma3)*tkm - 2.*k*(gamma3 -
                                                                alpha2*cosine_zenith)*tm)  # Equation 14.

    T = tm*(1. - (single_scatter_albedo/((1. - k*k*cosine_zenith*cosine_zenith) *
                                         ((k + gamma1)*tkp + (k - gamma1)*tkm))) *
            ((1. + k*cosine_zenith) *
             (alpha1 + k*gamma4)*tkp - (1. - k*cosine_zenith)*(alpha1 - k*gamma4)*tkm -
             2.*k*(gamma4 + alpha1*cosine_zenith)*tp))  # Equation 15.
    return R, T, T_pure


def scaling(optical_depth, single_scatter_albedo, asymmetry_factor):
    """Performs the scaling required by the delta-Eddington method.  The scaling
       parameters are described in https://doi.org/10.1175/1520-0469(1976)033<2452:TDEAFR>2.0.CO;2.

    Args:
        optical_depth: Optical depth.
        single_scatter_albedo: Single-scatter albedo.
        asymmetry_factor: Asymmetry factor.

    Returns:
        tau: Scaled optical depth.
        omega: Scaled single-scatter albedo.
        g: Scaled asymmetry factor.
    """
    g = asymmetry_factor/(asymmetry_factor + 1.)  # Equation 5b.
    f = g*g  # Equation 5a.
    omega = (1. - f)*single_scatter_albedo/(1. - single_scatter_albedo*f)  # Equation 14.
    tau = optical_depth*(1. - single_scatter_albedo*f)  # Equation 13.
    return tau, omega, g
