using Optim, Roots



function kohler_curve(d_dry,
    κ,
    σ = 0.072225,
    d_wet = 1,
    ρ_w = 997,
    M_w = 0.018015,
    R = 8.314,
    T = 298
    )

    """
    This function calculates a kohler curve based on the
    parameterization given by Petters & Kreidenweiss (2007)

    Non-optional arguments:
        d_dry [nm]
        κ, hygroscopicty parameter

    Optional arguments:
        σ, surface tension of water       [N m-1]
        ρ_w, density of water             [kg m-2]
        M_w, mmolecular weight of water   [kg m-3]
        R, gas constant                   [J mol-1 K-1]
        T, temperature                    [K]
        d_wet, wet diameters              [nm]
    """

    # if you don't specify the wet diameters then the curve is
    # generated from the dry diameter to 5 μm
    if d_wet == 1
        d_wet = collect(d_dry:10000)
    end

    # Calculate the supersaturation profile
    # Note on units: d_dry and d_wet are both presented in [nm].
    # The water activity term takes these like units, but the
    # exponential argument takes SI units, so convert d_wet to [m].
    S = @. (d_wet^3 - d_dry^3) / (d_wet^3 - d_dry^3 * (1 - κ)) * exp(4 * σ * M_w / (R * T * ρ_w * d_wet/10^9))

    return S
end



function calc_kappa(S_crit, d_dry)

    """
    This is a non-elegant function to estimate the
    hygroscypocity parameter, κ, base do on the parameterization
    given in Petters & Kreidenweiss (2007). This is intended
    for using the laboratory results. Requires the kohler curve
    calculation function.

    Arguments (both required):
        d_dry = dry diameter before hygroscopic growth [nm]
        S_crit = critical supersaturation; aka where AF = 0.5
    """

    # We will assume the corresponding critical diameter is < 5 μm and
    # calculate for wet diameters less than that cutoff.
    d_wet = collect(d_dry:0.1:5000);

    # The input we use is the critical supersaturation as a percent.
    # let's convert this to the useful supersation ratio
    S_crit = S_crit/100 + 1;

    # Basically, there is a unique critical supersaturation and critical
    # diameter for every value of kappa, so we're just gonna grid search for it.
    # It's possible kappa falls out of this range. 4 sig figs is probs fine
    κ_range = collect(0.01:0.00001:1.4);

    # To determine when the grid search ends, we will compare adjacent values.
    # changing kappa changes S_crit monotonically in one direction. We will create
    # these variables to help us store the previous values
    S_previous = 0;

    for i = 1:length(κ_range)

        # calculate the kohler curve
        S = kohler_curve(d_dry, κ_range[i]);

        # extract the calculated saturation ratio
        estimated_S_crit = maximum(S);

        # The difference betyween the true value, S_crit, and that we have
        # determined from kappa should decrease until we reach the max. At
        # that point it will begin to increase. Right then, we have found κ
        if abs((S_crit - estimated_S_crit)/S_crit) > abs((S_crit - S_previous)/S_crit)
            global κ_final           # allows us to pull κ from the loop
            κ_final = κ_range[i]
            break
        else
            S_previous = estimated_S_crit
        end
    end

    return κ_final
end



function calc_A_aw(d_wet,
        d_seed,
        d_coat,
        A0,
        C0,
        m_σ,
        ρ_inorg,
        M_inorg,
        VHF_inorg,
        ρ_org,
        M_org,
        VHF_org,
        ρ_w = 997,
        M_w = 0.018015,
        N_A = 6.022e23,
        R = 8.314,
        T = 298
    )

    """
    This function estimates A, used for calculating σ,
    given by the compressed film model. Further it estiates
    the water acitvity directly without the need to use κ.
    It uses a root finder to solve equations (8) - (10)
    from Forestieri et al. (2018). I went through and
    tried to solve somewhat analytically, hence the constants
    as defined here. Ultimately, I just tried to simplify the
    equations for the bisection root finding.

    Package dependence: Rooots.jl
    """

    # quick caluclation of molar volume [m3 mol-1]
    v_org = M_org / ρ_org;
    v_w = M_w / ρ_w;

    # Defining the constants
    k1 = m_σ * N_A / (2 * R * T);                                            # from equation (8)
    k2 = (d_coat^3 - d_seed^3) / (d_wet^3 * v_org);                          # from equation (9)
    k3 = (6 * v_org * (d_wet)^2) / (N_A * (d_coat^3 - d_seed^3) * 10^-9);    # from equation (10), note unit conversion for [nm] --> [m]
#     @show k1, k2, k3, A0, C0

    # Define function on which to find the root.
    # Here, x = f_surf
    f(x) = C0*exp((A0^2 - k3^2 / x^2) * k1 ) - (1-x)*k2;
#     D(f) = x -> ForwardDiff.derivative(f,float(x))
#     x0 = 0.5
#     f_surf = find_zero((f,D(f)), x0, Roots.Newton())
    # Solving for the root. f_surf has to be between 0 and 1
    # I think, so hopefully there isn't multiple roots lol
    f_surf = find_zero(f, (-0.001, 1.001))#, Bisection());
#       @show f_surf

    # Now can calculate A
    A = k3 / f_surf;

    # Finally, estimate the water activity as the fractional
    # molar amount of water in the solution, which has organic,
    # inorganic, and water as its constituents. I don't convert
    # the diameter units becasue it should cancel in the final
    # estimation of water activity.
    mol_water = pi/6 * ((d_wet)^3 -(d_coat)^3) * ρ_w / M_w;
    mol_inorg = pi/6 * (d_seed)^3 * ρ_inorg / M_inorg;
    mol_org = (1 - f_surf) * pi/6 * ((d_coat)^3 -(d_seed)^3) * ρ_org / M_org;
    a_w = mol_water / (mol_water + VHF_inorg*mol_inorg + VHF_org*mol_org);
#     @show d_wet, d_seed, d_coat#, ρ_w / M_w, ρ_inorg / M_inorg, ρ_org / M_org
#     @show mol_water, mol_inorg, mol_org
    return A, a_w, f_surf
end



function calc_sigma(A,
        A0,
        m_σ,
        σ_min = 0.04,
        σ_water = 0.072225,
    )

    """
    This function estimates σ given by the
    compressed film model. This is given by equation (7)
    in Forestieri et al. (2018).
    """

    σ_obs = min(σ_water, max(σ_water - (A0 - A)*m_σ, σ_min));

    return σ_obs
end



function calc_saturation_ratio(a_w,
        σ,
        d_wet,
        M_w = 0.018015,
        R = 8.314,
        T = 298,
        ρ_w = 997
    )

    """
    This function estimates the saturation ratio
    based on the Kohler Equation for a droplet
    containing solutes, Eq. (17.10) in Seinfeld &
    Pandis.

    Input parameters
        a_w, calculated water activity
        σ, dynamic surface tension [N/m]
        d_wet, wet diameter [nm]
    """

    kelvin_term = exp(4 * σ * M_w / (R * T * ρ_w * d_wet/10^9));
    S = @. a_w * kelvin_term;

    return S
end
