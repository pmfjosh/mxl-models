import numpy as np
from mxlpy import Derived, InitialAssignment, Model, Variable


def value(x: float) -> float:
    return x


def one_div(x: float) -> float:
    return 1.0 / x


def moiety_1(concentration: float, total: float) -> float:
    return total - concentration


def div(x: float, y: float) -> float:
    return x / y


def co2_initial(ca, Kh_co2):
    return 0.3 * ca / Kh_co2


def ci_initial(ca):
    return 0.65 * ca


def _pi_bellasio2019(total, pga, dhap, ru5p, rubp, atp):
    return total - pga - dhap - ru5p - 2 * rubp - atp


def _Et(vmax_rub, kcat_rub, V_m):
    return (vmax_rub / kcat_rub) / V_m


def _km_rubp_extra(pga, nadp, adp, pi, km_rubp, ki_pga, ki_nadp, ki_adp, ki_pi):
    return km_rubp * (1 + pga / ki_pga + nadp / ki_nadp + adp / ki_adp + pi / ki_pi)


def _f_rubp(rubp, Et, k_extra_rubp):
    top = (
        Et
        + k_extra_rubp
        + rubp
        - np.sqrt((Et + k_extra_rubp + rubp) ** 2 - 4 * rubp * Et)
    )
    bottom = 2 * Et
    return top / bottom


def non_rect_hyperbole(x, alpha, V0, theta):
    # print(np.sqrt((alpha * x + 1 - V0)**2 - 4 * alpha * x * theta))
    # top = alpha * x + 1 - V0 - np.sqrt((alpha * x + 1 - V0)**2 - 4 * alpha * x * theta)
    # bottom = 2 * theta
    return (
        (alpha * x + 1 - V0) / (2 * theta)
        - np.sqrt((alpha * x + 1 - V0) ** 2 - 4 * alpha * x * theta * (1 - V0))
        / (2 * theta)
        + V0
    )


def _Ract_eq(co2, ppfd, alpha_ppfd, V0_ppfd, theta_ppfd, alpha_co2, V0_co2, theta_co2):
    f_ppfd = non_rect_hyperbole(ppfd, alpha_ppfd, V0_ppfd, theta_ppfd)
    f_co2 = non_rect_hyperbole(co2, alpha_co2, V0_co2, theta_co2)
    return f_ppfd * f_co2


def _i20(pfd, s):
    return pfd * s


def _i10(i_20, y_ii_ll, y_i_ll):
    return i_20 * y_ii_ll / y_i_ll


def _chi(f_cyc, y_ii_ll):
    return f_cyc / (1 + f_cyc + y_ii_ll)


def _i1(chi, i10):
    return (1 + chi) * i10


def _f_cyc(j_atp, j_nadph, v_atp, v_fnr):
    return max(0, -1 + 15 ** (v_atp / j_atp - v_fnr / j_nadph))


def _i2(y_ii_ll, chi, i20):
    return (1 / y_ii_ll - chi) * i20 * y_ii_ll


def _y_ii(y_ii_ll, v_atp, j_atp, v_fnr, j_nadph, pfd, alpha, V0, theta):
    f_ppfd = non_rect_hyperbole(pfd, alpha, V0, theta)
    return y_ii_ll * (v_atp / j_atp) * (v_fnr / j_nadph) * (1 - max(0, f_ppfd))


def _j2(i2, y_ii):
    return i2 * y_ii


def _j1(j2, f_cyc):
    return j2 / 1 - f_cyc


def _f_pseudocyc(j_nadph, o2, v_fnr, f_pseudocycNR):
    return f_pseudocycNR + 4 * o2 * (1 - v_fnr / j_nadph)


def _j_nadph_steady(j1, f_cyc, f_pseudocyc):
    top = 1 - f_cyc - f_pseudocyc
    bottom = 2
    return (j1 * top / bottom) / 1000  # Added from Excel


def _j_atp_steady(j2, j1, f_cyc, fq, f_ndh, h):
    jcyt = (1 - fq) * j1
    jq = fq * j1
    jndh = f_cyc * f_ndh * j1

    return ((j2 + jcyt + 2 * jq + 2 * jndh) / h) / 1000  # Added from Excel


def _gs_steady(tau0, f_rubp, chi_beta, phi, pi_e, Kh, Ds, gs0):

    tau = tau0 + f_rubp
    top = chi_beta * tau * (phi + pi_e)
    bottom = 1 + chi_beta * tau * (1 / Kh) * Ds

    return max(gs0, top / bottom)


def _calc_ass(vc, vo, RLight):
    return vc - 0.5 * vo - RLight


def include_derived_quantities(m: Model):

    m.add_derived(
        name="ADP_st",
        fn=moiety_1,
        args=["ATP_st", "AP_tot"],
    )

    m.add_derived(
        name="NADP_st",
        fn=moiety_1,
        args=["NADPH_st", "NADP_tot"],
    )

    m.add_derived(
        name="Pi_st",
        fn=_pi_bellasio2019,
        args=["Pi_tot", "PGA", "DHAP", "RU5P", "RUBP", "ATP_st"],
    )

    m.add_derived(
        name="Et",
        fn=_Et,
        args=["vmax_v_RuBisCO_c", "kcat_v_RuBisCO_c", "V_m"],
    )

    m.add_derived(
        name="km_v_RuBisCO_c_RUBP_extra",
        fn=_km_rubp_extra,
        args=[
            "PGA",
            "NADP_st",
            "ADP_st",
            "Pi_st",
            "km_v_RuBisCO_c_RUBP",
            "ki_v_RuBisCO_c_PGA",
            "ki_v_RuBisCO_c_NADP_st",
            "ki_v_RuBisCO_c_ADP_st",
            "ki_v_RuBisCO_c_Pi_st",
        ],
    )

    m.add_derived(
        name="f_rubp",
        fn=_f_rubp,
        args=["RUBP", "Et", "km_v_RuBisCO_c_RUBP_extra"],
    )

    m.add_derived(
        name="O2",
        fn=div,
        args=["p_o2", "Kh_o2"],
    )

    m.add_derived(
        name="Ract_eq",
        fn=_Ract_eq,
        args=[
            "CO2",
            "PPFD",
            "alpha_ppfd_rub",
            "V0_ppfd_rub",
            "theta_ppfd_rub",
            "alpha_co2",
            "V0_co2",
            "theta_co2",
        ],
    )

    m.add_derived(
        name="I2_0",
        fn=_i20,
        args=["PPFD", "s"],
    )

    m.add_derived(
        name="I1_0",
        fn=_i10,
        args=["I2_0", "PhiPSII_LL", "PhiPSI_LL"],
    )

    m.add_derived(
        name="chi",
        fn=_chi,
        args=["f_cyc", "PhiPSII_LL"],
    )

    m.add_derived(
        name="I1",
        fn=_i1,
        args=["chi", "I1_0"],
    )

    m.add_derived(
        name="f_cyc",
        fn=_f_cyc,
        args=["J_ATP", "J_NADPH", "v_ATPsynth", "v_FNR"],
    )

    m.add_derived(
        name="I2",
        fn=_i2,
        args=["PhiPSII_LL", "chi", "I2_0"],
    )

    m.add_derived(
        name="PhiPSII",
        fn=_y_ii,
        args=[
            "PhiPSII_LL",
            "v_ATPsynth",
            "J_ATP",
            "v_FNR",
            "J_NADPH",
            "PPFD",
            "alpha_ppfd_PhiPSII",
            "V0_ppfd_PhiPSII",
            "theta_ppfd_PhiPSII",
        ],
    )

    m.add_derived(
        name="J2",
        fn=_j2,
        args=["I2", "PhiPSII"],
    )

    m.add_derived(
        name="J1",
        fn=_j1,
        args=["J2", "f_cyc"],
    )

    m.add_derived(
        name="f_pseudocyc",
        fn=_f_pseudocyc,
        args=["J_NADPH", "O2", "v_FNR", "f_pseudocycNR"],
    )

    m.add_derived(
        name="J_NADPH_steady",
        fn=_j_nadph_steady,
        args=["J1", "f_cyc", "f_pseudocyc"],
    )

    m.add_derived(
        name="J_ATP_steady",
        fn=_j_atp_steady,
        args=["J2", "J1", "f_cyc", "fq", "f_ndh", "h"],
    )

    m.add_derived(
        name="gs_steady",
        fn=_gs_steady,
        args=["tau0", "f_rubp", "chi_beta", "phi", "pi_e", "Kh", "Ds", "gs0"],
    )

    m.add_derived(
        name="A",
        fn=_calc_ass,
        args=["v_RuBisCO_c", "rubisco_oxygenase", "RLight"],
    )

    return m


def ract_gs_time_dependance(x, x_steady, inc, dec):
    if x < x_steady:
        return (x_steady - x) / inc
    else:
        return (x_steady - x) / dec


def atp_nadph_time_dependance(j_x, j_x_steady, kj_x):
    if j_x < j_x_steady:
        return (j_x_steady - j_x) / kj_x
    else:
        return (j_x_steady - j_x) / (0.1 * kj_x)


def _rubisco_carboxylation_bellasio(
    rubp, co2, Ract, km_co2, o2, km_o2, vmax_rc, f_rubp, k_extra_rubp
):
    k_extra_co2 = km_co2 * (1 + o2 / km_o2)

    top = vmax_rc * Ract * f_rubp * rubp * co2
    bottom = (k_extra_co2 + co2) * (k_extra_rubp + rubp)

    return top / bottom


def neg_one_div(x: float) -> float:
    return -1.0 / x


def two_div(x: float) -> float:
    return 2.0 / x


def _rubisco_oxygenase_bellasio(co2, o2, S_co_gas, v_c, Kh_o2, Kh_co2):
    S_co_liq = S_co_gas / Kh_o2 * Kh_co2
    gamma_star = 1 / (2 * S_co_liq)

    return v_c * 2 * gamma_star * o2 / co2


def neg_half_div(x: float) -> float:
    return -0.5 / x


def half_div(x: float) -> float:
    return 0.5 / x


def _prkase(
    atp,
    rubp,
    ru5p,
    pga,
    adp,
    pi,
    vmax,
    k_eq,
    km_atp,
    ki_adp,
    km_ru5p,
    ki_pga,
    ki_rubp,
    ki_pi,
):
    top = vmax * (atp * ru5p - adp * rubp / k_eq)  # ERROR: Different from paper (typo?)
    bottom = (atp + km_atp * (1 + adp / ki_adp)) * (
        ru5p + km_ru5p * (1 + pga / ki_pga + rubp / ki_rubp + pi / ki_pi)
    )
    return top / bottom


def neg_fivethirds_div(x):
    return -(5 / 3) * (1 / x)


def _v_pgareduction(atp, pga, nadph, adp, vmax, km_atp, km_pga, km_nadph, ki_adp):
    top = vmax * atp * pga * nadph
    bottom = (
        (pga + km_pga * (1 + adp / ki_adp))
        * (atp + km_atp * (1 + adp / ki_adp))
        * (nadph + km_nadph * (1 + adp / ki_adp))
    )
    return top / bottom


def _v_carbohydrate_synthesis(
    dhap, pi, adp, vmax, v_pgareduction, keq, km_dhap, ki_adp
):
    top = vmax * (dhap - 0.4) * (1 - np.abs(v_pgareduction) * pi / keq)
    bottom = dhap + km_dhap * (1 + adp / ki_adp)
    return top / bottom


def _v_rpp(dhap, ru5p, vmax, k_eq, km_dhap):
    top = vmax * (dhap - ru5p / k_eq)
    bottom = dhap + km_dhap
    return top / bottom


def _v_co2_hydration(co2, hco3, proton, vmax, k_eq, km_co2, km_hco3):
    top = vmax * (co2 - hco3 * proton / k_eq)
    bottom = km_co2 * (1 + co2 / km_co2 + hco3 / km_hco3)
    return top / bottom


def neg_onethirds_div(x):
    return -(1 / 3) * (1 / x)


def _v_fnr(nadph, nadp, j_nadph, k_eq, km_nadp, km_nadph):
    top = j_nadph * (nadp - nadph / k_eq)
    bottom = km_nadp * (1 + nadp / km_nadp + nadph / km_nadph)
    return top / bottom


def _v_atp(atp, adp, pi, j_atp, k_eq, km_adp, km_pi, km_atp):
    top = j_atp * (adp * pi - atp / k_eq)
    bottom = (
        km_adp
        * km_pi
        * (1 + adp / km_adp + atp / km_atp + pi / km_pi + adp * pi / (km_adp * km_pi))
    )
    return top / bottom


def _co2_diss(ci, co2, gm, Kh_co2):
    return (gm * (ci - co2 * Kh_co2)) / 1000


def _stom_diff(ci, gs, ca):
    return (gs * (ca - ci)) / 1000


def include_rates(m: Model):

    m.add_reaction(
        name="Ract_rate",
        fn=ract_gs_time_dependance,
        args=["Ract", "Ract_eq", "tau_i", "tau_d"],
        stoichiometry={
            "Ract": 1,
        },
    )
    m.add_reaction(
        name="v_J_NADPH",
        fn=atp_nadph_time_dependance,
        args=["J_NADPH", "J_NADPH_steady", "Kj_NADPH"],
        stoichiometry={
            "J_NADPH": 1,
        },
    )
    m.add_reaction(
        name="v_J_ATP",
        fn=atp_nadph_time_dependance,
        args=["J_ATP", "J_ATP_steady", "Kj_ATP"],
        stoichiometry={
            "J_ATP": 1,
        },
    )
    m.add_reaction(
        name="v_gs",
        fn=ract_gs_time_dependance,
        args=["gs", "gs_steady", "Ki", "Kd"],
        stoichiometry={
            "gs": 1,
        },
    )
    m.add_reaction(
        name="v_RuBisCO_c",
        fn=_rubisco_carboxylation_bellasio,
        args=[
            "RUBP",
            "CO2",
            "Ract",
            "km_v_RuBisCO_c_CO2",
            "O2",
            "km_v_RuBisCO_c_O2",
            "vmax_v_RuBisCO_c",
            "f_rubp",
            "km_v_RuBisCO_c_RUBP_extra",
        ],
        stoichiometry={
            "CO2": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "RUBP": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "PGA": Derived(fn=two_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="rubisco_oxygenase",
        fn=_rubisco_oxygenase_bellasio,
        args=["CO2", "O2", "S_co_gas", "v_RuBisCO_c", "Kh_o2", "Kh_co2"],
        stoichiometry={
            "RUBP": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "PGA": Derived(fn=one_div, args=["V_m"], unit=None),
            "ATP_st": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "NADPH_st": Derived(fn=neg_half_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="glycine_decarboxylase",
        fn=value,
        args=["rubisco_oxygenase"],
        stoichiometry={
            "CO2": Derived(fn=half_div, args=["V_m"], unit=None),
            "PGA": Derived(fn=half_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_PRKase",
        fn=_prkase,
        args=[
            "ATP_st",
            "RUBP",
            "RU5P",
            "PGA",
            "ADP_st",
            "Pi_st",
            "vmax_v_PRKase",
            "keq_v_PRKase",
            "km_v_PRKase_ATP_st",
            "ki_v_PRKase_ADP_st",
            "km_v_PRKase_RU5P",
            "ki_v_PRKase_PGA",
            "ki_v_PRKase_RUBP",
            "ki_v_PRKase_Pi_st",
        ],
        stoichiometry={
            "RUBP": Derived(fn=one_div, args=["V_m"], unit=None),
            "DHAP": Derived(fn=neg_fivethirds_div, args=["V_m"], unit=None),
            "ATP_st": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "RU5P": Derived(fn=neg_one_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_pgareduction",
        fn=_v_pgareduction,
        args=[
            "ATP_st",
            "PGA",
            "NADPH_st",
            "ADP_st",
            "vmax_v_pgareduction",
            "km_v_pgareduction_ATP_st",
            "km_v_pgareduction_PGA",
            "km_v_pgareduction_NADPH_st",
            "ki_v_pgareduction_ADP_st",
        ],
        stoichiometry={
            "PGA": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "DHAP": Derived(fn=one_div, args=["V_m"], unit=None),
            "ATP_st": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "NADPH_st": Derived(fn=neg_one_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_carbohydrate_synthesis",
        fn=_v_carbohydrate_synthesis,
        args=[
            "DHAP",
            "Pi_st",
            "ADP_st",
            "vmax_v_carbohydrate_synthesis",
            "v_pgareduction",
            "keq_v_carbohydrate_synthesis",
            "km_v_carbohydrate_synthesis_DHAP",
            "ki_v_carbohydrate_synthesis_ADP_st",
        ],
        stoichiometry={
            "DHAP": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "ATP_st": Derived(fn=neg_half_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_rpp",
        fn=_v_rpp,
        args=["DHAP", "RU5P", "vmax_v_rpp", "keq_v_rpp", "km_v_rpp_DHAP"],
        stoichiometry={
            "RU5P": Derived(fn=one_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_co2_hydration",
        fn=_v_co2_hydration,
        args=[
            "CO2",
            "HCO3",
            "H",
            "vmax_v_co2_hydration",
            "keq_v_co2_hydration",
            "km_v_co2_hydration_CO2",
            "km_v_co2_hydration_HCO3",
        ],
        stoichiometry={
            "CO2": Derived(fn=neg_one_div, args=["V_m"], unit=None),
            "HCO3": Derived(fn=one_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_RLight",
        fn=value,
        args=["RLight"],
        stoichiometry={
            "CO2": Derived(fn=one_div, args=["V_m"], unit=None),
            "PGA": Derived(fn=neg_onethirds_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_FNR",
        fn=_v_fnr,
        args=[
            "NADPH_st",
            "NADP_st",
            "J_NADPH",
            "keq_v_FNR",
            "km_v_FNR_NADP_st",
            "km_v_FNR_NADPH_st",
        ],
        stoichiometry={
            "NADPH_st": Derived(fn=one_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="v_ATPsynth",
        fn=_v_atp,
        args=[
            "ATP_st",
            "ADP_st",
            "Pi_st",
            "J_ATP",
            "keq_v_ATPsynth",
            "km_v_ATPsynth_ADP_st",
            "km_v_ATPsynth_Pi_st",
            "km_v_ATPsynth_ATP_st",
        ],
        stoichiometry={
            "ATP_st": Derived(fn=one_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="CO2_dissolution",
        fn=_co2_diss,
        args=["Ci", "CO2", "gm", "Kh_co2"],
        stoichiometry={
            "Ci": -1,
            "CO2": Derived(fn=one_div, args=["V_m"], unit=None),
        },
    )
    m.add_reaction(
        name="CO2_stomatal_diffusion",
        fn=_stom_diff,
        args=["Ci", "gs", "Ca"],
        stoichiometry={
            "Ci": 1,
        },
    )

    return m


def get_bellasio_2019() -> Model:
    m = Model()

    m.add_parameters(
        {
            "AP_tot": 1.5,
            "Pi_tot": 15,
            "p_o2": 210000,
            "Kh_o2": 833300,
            "V_m": 0.03,
            "PPFD": 1500,
            "RLight": 0.001,
            "s": 0.43,
            "PhiPSII_LL": 0.72,
            "PhiPSI_LL": 1,
            "alpha_ppfd_PhiPSII": 0.00125,
            "V0_ppfd_PhiPSII": -0.8,
            "theta_ppfd_PhiPSII": 0.7,
            "f_pseudocycNR": 0.01,
            "fq": 1,
            "f_ndh": 0,
            "h": 4,
            "Ca": 350,
            "alpha_ppfd_rub": 0.0018,
            "V0_ppfd_rub": 0.16,
            "theta_ppfd_rub": 0.95,
            "alpha_co2": 400,
            "V0_co2": -0.2,
            "theta_co2": 0.98,
            "tau_i": 360,
            "tau_d": 1200,
            "km_v_RuBisCO_c_CO2": 0.014,
            "km_v_RuBisCO_c_RUBP": 0.02,
            "km_v_RuBisCO_c_O2": 0.222,
            "ki_v_RuBisCO_c_PGA": 2.52,
            "ki_v_RuBisCO_c_NADP_st": 0.21,
            "ki_v_RuBisCO_c_ADP_st": 0.2,
            "ki_v_RuBisCO_c_Pi_st": 3.6,
            "vmax_v_RuBisCO_c": 0.2,
            "kcat_v_RuBisCO_c": 4.7,
            "S_co_gas": 2200,
            "vmax_v_PRKase": 1.17,
            "keq_v_PRKase": 6846,
            "km_v_PRKase_ATP_st": 0.625,
            "ki_v_PRKase_ADP_st": 0.1,
            "km_v_PRKase_RU5P": 0.034,
            "ki_v_PRKase_PGA": 2,
            "ki_v_PRKase_RUBP": 0.7,
            "ki_v_PRKase_Pi_st": 4,
            "vmax_v_pgareduction": 0.4,
            "km_v_pgareduction_ATP_st": 0.3,
            "km_v_pgareduction_PGA": 10,
            "km_v_pgareduction_NADPH_st": 0.05,
            "ki_v_pgareduction_ADP_st": 0.89,
            "vmax_v_carbohydrate_synthesis": 0.2235,
            "keq_v_carbohydrate_synthesis": 0.8,
            "km_v_carbohydrate_synthesis_DHAP": 22,
            "ki_v_carbohydrate_synthesis_ADP_st": 1,
            "vmax_v_rpp": 0.0585,
            "keq_v_rpp": 0.06,
            "km_v_rpp_DHAP": 0.5,
            "H": 5.012e-05,
            "vmax_v_co2_hydration": 200,
            "keq_v_co2_hydration": 0.00056,
            "km_v_co2_hydration_CO2": 2.8,
            "km_v_co2_hydration_HCO3": 34,
            "keq_v_FNR": 502,
            "km_v_FNR_NADP_st": 0.0072,
            "km_v_FNR_NADPH_st": 0.036,
            "Kj_NADPH": 200,
            "keq_v_ATPsynth": 5734,
            "km_v_ATPsynth_ADP_st": 0.014,
            "km_v_ATPsynth_Pi_st": 0.3,
            "km_v_ATPsynth_ATP_st": 0.3,
            "Kj_ATP": 200,
            "gm": 0.5,
            "Kh_co2": 30303.0303030303,
            "Kd": 150,
            "Ki": 900,
            "tau0": -0.1,
            "chi_beta": 0.5,
            "phi": 0,
            "pi_e": 1.2,
            "Kh": 12,
            "Ds": 10,
            "gs0": 0.01,
            "NADP_tot": 0.5,
        }
    )

    m.add_variables(
        {
            "CO2": InitialAssignment(fn=co2_initial, args=["Ca", "Kh_co2"]),
            "HCO3": Variable(0.1327),
            "RUBP": Variable(2),
            "PGA": Variable(4),
            "DHAP": Variable(4),
            "ATP_st": Variable(0.68),
            "NADPH_st": Variable(0.21),
            "RU5P": Variable(0.34),
            "Ract": Variable(1),
            "J_NADPH": Variable(0.1),
            "J_ATP": Variable(0.16),
            "Ci": InitialAssignment(fn=ci_initial, args=["Ca"]),
            "gs": Variable(0.334934046786077),
        }
    )

    m = include_derived_quantities(m)
    m = include_rates(m)

    return m
