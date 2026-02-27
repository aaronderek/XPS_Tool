"""
M1: Triangular potential well model (constant field approximation).

Physical picture: Constant electric field at the surface.
Potential: V(z) = q·Es·z for 0 ≤ z ≤ W
           V(z) = 0 for z > W
"""

import numpy as np
from physics.constants import Q, EPS0, HBAR, M2_TO_CM2
from physics.units import nm_to_m, eV_to_J, J_to_eV, ns_to_display, ns_from_display


class TriangularModel:
    """M1: Triangular potential well with constant electric field."""

    def __init__(self, m_star, epsilon_r, W_nm):
        """
        Initialize the triangular model.

        Parameters
        ----------
        m_star : float
            Effective mass in kg
        epsilon_r : float
            Relative permittivity
        W_nm : float
            Depletion width in nanometers
        """
        self.m_star = m_star
        self.epsilon_r = epsilon_r
        self.W_nm = W_nm
        self.W = nm_to_m(W_nm)  # Convert to meters
        self.epsilon = EPS0 * epsilon_r

    def calculate_ns(self, Phi_s_eV):
        """
        Calculate sheet density from surface potential.

        Parameters
        ----------
        Phi_s_eV : float or array
            Surface potential in eV

        Returns
        -------
        ns : float or array
            Sheet density in m⁻²
        """
        # Note: Phi_s_eV is the surface potential. For electrons:
        # "1 eV potential" means electric potential of 1 V
        # Phi_s = Es * W  =>  Es = Phi_s / W
        Phi_s_V = Phi_s_eV  # Numeric value is the same
        Es = Phi_s_V / self.W

        # Gauss's law: ns = (epsilon * Es) / q
        ns = (self.epsilon * Es) / Q

        return ns

    def calculate_Phi_s(self, ns):
        """
        Calculate surface potential from sheet density (inverse calculation).

        Parameters
        ----------
        ns : float or array
            Sheet density in m⁻²

        Returns
        -------
        Phi_s_eV : float or array
            Surface potential in eV
        """
        # ns = (epsilon * Es) / q  =>  Es = (ns * q) / epsilon
        Es = (ns * Q) / self.epsilon

        # Phi_s = Es * W
        Phi_s_V = Es * self.W

        return Phi_s_V  # Return as eV (numerically equals V)

    def calculate_Delta_WF(self, ns):
        """
        Calculate work function change from sheet density.

        Parameters
        ----------
        ns : float or array
            Sheet density in m⁻²

        Returns
        -------
        Delta_WF : float or array
            Work function change in eV (without adsorbates)
        """
        # Delta_WF = -Phi_s = -(q*W/epsilon) * ns
        Phi_s_eV = self.calculate_Phi_s(ns)
        return -Phi_s_eV

    def get_surface_field(self, Phi_s_eV):
        """
        Get surface electric field.

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV

        Returns
        -------
        Es : float
            Surface electric field in V/m
        """
        Phi_s_V = Phi_s_eV  # Potential in volts
        return Phi_s_V / self.W

    def get_potential(self, Phi_s_eV, z_array):
        """
        Get potential profile V(z).

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV
        z_array : array
            Array of z positions in meters

        Returns
        -------
        V_z : array
            Potential energy in Joules at each z position
        """
        Phi_s_V = Phi_s_eV  # Potential in volts

        # Electric potential: V_elec(z) = -Phi_s * (1 - z/W) for z ≤ W
        # (Negative because we're measuring from bulk, surface is lower)
        z_ratio = z_array / self.W
        V_elec = -Phi_s_V * (1 - z_ratio)

        # Potential energy for electron: U = -e * V_elec = +e * Phi_s * (1-z/W)
        # (Positive because electron is negatively charged)
        V_z = -Q * V_elec  # This gives positive energy (barrier for electrons)

        # Set potential to 0 beyond W
        V_z[z_array > self.W] = 0

        return V_z

    def get_subband_energies(self, Phi_s_eV, n_levels=3):
        """
        Calculate quantum subband energy levels.

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV
        n_levels : int
            Number of subband levels to calculate

        Returns
        -------
        energies : array
            Subband energies in eV
        """
        Es = self.get_surface_field(Phi_s_eV)

        # Airy function solutions: En = an * (ℏ²/(2m*))^(1/3) * (q*Es)^(2/3)
        # Numerical constants for first few levels
        a_n = [2.338, 4.088, 5.521, 6.787, 7.944]

        energies = []
        for n in range(min(n_levels, len(a_n))):
            En = a_n[n] * (HBAR**2 / (2 * self.m_star))**(1/3) * (Q * Es)**(2/3)
            energies.append(J_to_eV(En))

        return np.array(energies)

    def get_slope(self):
        """
        Get the slope of Delta_WF vs ns curve.

        Returns
        -------
        slope : float
            Slope in eV/(10¹³ cm⁻²)
        """
        # Slope = -(q*W/epsilon) in SI units
        # slope_SI has units V·m² (potential times area)
        slope_V_m2 = -(Q * self.W) / self.epsilon  # V·m²

        # Convert to eV/(10¹³ cm⁻²)
        # Unit analysis: 1 cm⁻² = (10⁻² m)⁻² = 10⁴ m⁻²
        # So: 10¹³ cm⁻² = 10¹³ × 10⁴ m⁻² = 10¹⁷ m⁻²
        # Therefore: ns[m⁻²] = ns_display[10¹³ cm⁻²] × 10¹⁷
        # From ΔWF = slope_V_m2 × ns = slope_display × ns_display:
        # slope_display = slope_V_m2 × (ns/ns_display) = slope_V_m2 × 10¹⁷
        slope_display = slope_V_m2 * 1e17  # eV/(10¹³ cm⁻²)

        return slope_display

    def __str__(self):
        return f"M1-Triangular(εᵣ={self.epsilon_r}, W={self.W_nm:.1f}nm)"
