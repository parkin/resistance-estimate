#!/usr/bin/env python
from scipy.optimize import fsolve
from scipy.integrate import quad
import math
import numpy as np
import matplotlib.pyplot as pp


def series(*resistances):
    """
    Returns the series resistance.
    """
    r_out = 0
    return sum(resistances)


def parallel(*resistances):
    """
    Returns the parallel resistances.
    """
    inverse = 0
    for resistance in resistances:
        inverse += 1.0/resistance

    return 1.0/inverse


def calculate_dv(r, v_macro, molecule_diameter, pore_height,
        pore_diameter, concentration_ratio):

    molecule_cross_section = math.pi * (molecule_diameter/2.0)**2
    numerator = 2.0 * v_macro * molecule_cross_section \
            * (4.0 * pore_height + pore_diameter) \
            * (concentration_ratio - 1.0)

    denominator = math.pi * math.log(concentration_ratio) \
            * (2.0 * pore_height + pore_diameter) \
            * (pore_diameter**2 * (concentration_ratio - 1)\
            + 4.0 * (2.0 * pore_height + pore_diameter) * r)

    return numerator / denominator


def calculate_dr(pore_diameter, pore_height, concentration_ratio,
        v_macro, molecule_diameter, resistivity, x_lim, y_lim,
        ids, gate_slope):
    multiplier = resistivity * gate_slope / ids

    pore_radius = pore_diameter / 2.0

    if pore_radius > x_lim and pore_radius > y_lim:
        return 0.0

    def inner_integrand(x, y):
        r = math.sqrt(x**2 + y*2)
        return calculate_dv(r, v_macro, molecule_diameter, pore_height,
                pore_diameter, concentration_ratio)

    def outer_integrand(y, x_lim):

        inner_integral = quad(inner_integrand, pore_radius, x_lim,
                args=(y))

        return 1.0 / inner_integral[0]

    outer_integral = quad(outer_integrand, pore_radius, y_lim,
            args=(x_lim))
    print("outer_integral: {0}".format(outer_integral))

    return multiplier * (1.0 / outer_integral[0])


def calculate_active_radius(pore_diameter, pore_height,
        concentration_ratio, v_macro, molecule_diameter, resistivity,
        ids, gate_slope, threshold=0.01):
    """
    Returns distance in nm where the dv multiplier falls below threshold.
    """
    print "molecule_diameter: ", molecule_diameter
    if molecule_diameter <= 0.0:
        return 0.0

    def func(r):
        dv = calculate_dv(r, v_macro, molecule_diameter,
                                    pore_height, pore_diameter,
                                    concentration_ratio)
        result = dv / v_macro - threshold
        return result
    result = fsolve(func, 0)
    if len(result) < 1 or result[0] < 1:
        return 0.0
    
    return result[0]


def calculate_r_active(pore_diameter, pore_height, 
        concentration_ratio, v_macro, molecule_diameter, resistivity,
        x, y, ids, gate_slope):
    """
    Calculates the resistance of the active region.
    :param x: 1/2 width
    :param y: 1/2 height
    """
    resistance = resistivity

    if 2.*x <= pore_diameter and 2.*y <= pore_diameter:
        return resistance

    dr = calculate_dr(pore_diameter, pore_height, concentration_ratio, v_macro, molecule_diameter, resistivity, x, y, ids, gate_slope)

    print("dr: {0:.2E}".format(dr))

    return resistance + dr


def calculate_r_out(length, width, active_radius, resistivity):
    return 2.0 * resistivity * length / (width - 2.0 * active_radius)


def calculate_r_pre(length, width, active_radius, resistivity):
    return resistivity * (length - 2.0 * active_radius) / (4.0 * active_radius)


def calculate_resistance(length, width, pore_diameter,
        pore_height, v_macro, resistivity, r_contact, ids,
        gate_slope, concentration_ratio, molecule_diameter=0.0):
    """

    :param length: Length of the channel.
    :param width: Width of the channel.
    :param pore_diameter: Diameter of the pore, in nm.
    :param pore_height: Height of the pore, in nm.
    :param v_macro: Macro voltage between cis and trans chambers. in V
    :param resistivity: Resistivity of the material. in V/A
    :param r_contact: Contact resistance. in Ohms
    :param ids: The drain-source current at device operation. in A
    :param gate_slope: The current/voltage slope of the gating curve 
        at operation. in A/V.
    :param concentration_ratio: Ratio of chamber concentrations Ccis/Ctrans.
    :param molecule_diameter: Diameter of the moledule, in nm. Default is dsDNA (2nm)

    The device looks like this.

    * Contact is the contact resistance.
    * Rout is the parallel portion on the edges of the device
        that are unaffected by any potential change
    * Rpre are the series portions that are unaffected by 
        any potential change
    * Rchannel is the portion that is affected by the potential change.


    |                Contact               |
    ----------------------------------------
    |    Rout    |    Rpre    |    Rout    |
    |            |            |            |
    |            |            |            |
    |            |------------|            |
    |            |  Rchannel  |            |
    |            | (pore is   |            |
    |            |   here)    |            |
    |            |------------|            |
    |            |            |            |
    |            |            |            |
    |            |    Rpre    |            |
    ----------------------------------------
    |                Contact               |

    """
    print("\n----calculate_resistance")

    active_radius = calculate_active_radius(pore_diameter, pore_height,
            concentration_ratio, v_macro, molecule_diameter, resistivity,
            ids, gate_slope)
    print("Active radius: {0:.2E}".format(active_radius))

    if 2*active_radius > width:
        x = length
    else:
        x = active_radius
    if 2*active_radius > length:
        y = length
    else:
        y = active_radius

    r_active = calculate_r_active(pore_diameter, 
            pore_height, concentration_ratio, v_macro, molecule_diameter,
            resistivity, x, y, ids, gate_slope)
    print("r_active: {0:.2E}".format(r_active))

    # Combine the resistance parts
    if 2*active_radius < length and active_radius > 0.0:
        r_pre = calculate_r_pre(length, width, active_radius, resistivity)
        print("r_pre: {0:.2E}".format(r_pre))
        r_sensitive = series(r_pre, r_active, r_pre)
    else:
        r_sensitive = r_active

    if 2*active_radius < width:
        r_out = calculate_r_out(length, width, active_radius, resistivity)
        print("r_out: {0:.2E}".format(r_out))
        r_channel = parallel(r_out, r_sensitive, r_out)
    else:
        r_channel = r_sensitive

    print("r_channel: {0:.2E}".format(r_channel))

    r_total = series(r_contact, r_channel, r_contact)

    return r_total


def calculate_percent_change(length, width, pore_diameter,
        pore_height, v_macro, resistivity, r_contact, ids,
        gate_slope, concentration_ratio, molecule_diameter=2.0):

    r_on = calculate_resistance(length, width, pore_diameter,
            pore_height, v_macro, resistivity, r_contact, ids,
            gate_slope, concentration_ratio, molecule_diameter)
    r_normal = calculate_resistance(length, width, pore_diameter,
            pore_height, v_macro, resistivity, r_contact, ids,
            gate_slope, concentration_ratio)
    print("\n")
    print("r_on: {0:.2E}, r_normal: {1:.2E}".format(r_on, r_normal))

    return (r_on - r_normal) / r_normal 

def do():
    length = 500.
    width = 200.
    pore_diameter = 3.0
    pore_height = 50.
    v_macro = 500.e-3
    ids = 1.e-11
    vds = 300.e-3
    resistivity = vds / ids
    gate_slope = 1.e-10 / 120.e-3
    concentration_ratio = 100.
    molecule_diameter = 2.
    r_contact = 1.e6

    print("******************************")
    print("length: {0} nm".format(length))
    print("width: {0} nm".format(width))
    print("pore_diameter: {0} nm".format(pore_diameter))
    print("pore_height: {0} nm".format(pore_height))
    print("molecule_diameter: {0} nm".format(molecule_diameter))
    print("v_macro: {0} V".format(v_macro))
    print("ids: {0} A".format(ids))
    print("vds: {0} V".format(vds))
    print("resistivity: {0:.2E} Ohms/square".format(resistivity))
    print("r_contact: {0:.2E} Ohms".format(r_contact))
    print("gate_slope: {0:.2E} A/V".format(gate_slope))
    print("concentration_ratio: {0}".format(concentration_ratio))
    print("******************************")

    ratio = calculate_percent_change(length, width, pore_diameter,
            pore_height, v_macro, resistivity, r_contact, ids,
            gate_slope, concentration_ratio, molecule_diameter)

    print("\n------------")
    print("Ratio: {0:.2E}".format(ratio))

    # Plot dv

    x = np.linspace(0, 1000, 100)
    y = np.empty(100)
    for i in range(len(x)):
        y[i] = calculate_dv(x[i], v_macro, molecule_diameter, pore_height,
                pore_diameter, concentration_ratio)
    pp.plot(x, y)
    pp.show()

if __name__ == "__main__":
    do()
