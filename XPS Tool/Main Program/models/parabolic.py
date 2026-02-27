"""
M3: Parabolic/Harmonic potential well model.

Physical picture: Electric field linearly decays to zero.
Potential: V(z) = -Φs·(1 - z/W)² for 0 ≤ z ≤ W
           V(z) = 0 for z > W
"""

import numpy as np
from physics.constants import Q, EPS0, HBAR, M2_TO_CM2
from physics.units import nm_to_m, eV_to_J, J_to_eV, ns_to_display, ns_from_display


class ParabolicModel:
    """M3: Parabolic potential well with linearly decaying field."""

    def __init__(self, m_star, epsilon_r, W_nm):
        """
        Initialize the parabolic model.

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
        # Es = 2*Phi_s / W  (surface field)
        Phi_s_V = Phi_s_eV  # Potential in volts
        Es = 2 * Phi_s_V / self.W

        # Gauss's law: ns = (epsilon * Es) / q
        # Substituting: ns = 2*epsilon*Phi_s / (q*W)
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
        # ns = 2*epsilon*Phi_s / (q*W)  =>  Phi_s = ns*q*W / (2*epsilon)
        Phi_s_V = (ns * Q * self.W) / (2 * self.epsilon)

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
        # Delta_WF = -Phi_s = -(q*W/(2*epsilon)) * ns
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
        return 2 * Phi_s_V / self.W

    def get_field(self, Phi_s_eV, z_array):
        """
        Get electric field profile E(z).

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV
        z_array : array
            Array of z positions in meters

        Returns
        -------
        E_z : array
            Electric field in V/m at each z position
        """
        Es = self.get_surface_field(Phi_s_eV)

        # E(z) = (2*Phi_s/W) * (1 - z/W) = Es * (1 - z/W)
        E_z = Es * (1 - z_array / self.W)

        # Set field to 0 beyond W
        E_z[z_array > self.W] = 0

        return E_z

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

        # Electric potential: V_elec(z) = -Phi_s * (1 - z/W)² for z ≤ W
        # (Negative because we're measuring from bulk, surface is lower)
        z_ratio = z_array / self.W
        V_elec = -Phi_s_V * (1 - z_ratio)**2

        # Potential energy for electron: U = -e * V_elec = +e * Phi_s * (1-z/W)²
        # (Positive because electron is negatively charged)
        V_z = -Q * V_elec  # This gives positive energy (barrier for electrons)

        # Set potential to 0 beyond W
        V_z[z_array > self.W] = 0

        return V_z

    def get_harmonic_levels(self, Phi_s_eV, n_levels=3):
        """
        Calculate approximate quantum levels using harmonic approximation.

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV
        n_levels : int
            Number of levels to calculate

        Returns
        -------
        energies : array
            Approximate energy levels in eV
        """
        # Near the surface, approximate as harmonic oscillator
        # ω = √(2*q*Es/(m*W))
        Es = self.get_surface_field(Phi_s_eV)
        omega = np.sqrt(2 * Q * Es / (self.m_star * self.W))

        # En = ℏ*ω*(n + 1/2)
        energies = []
        for n in range(n_levels):
            En = HBAR * omega * (n + 0.5)
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
        # Slope = -(q*W/(2*epsilon)) in SI units
        # This is exactly half of M1's slope
        slope_V_m2 = -(Q * self.W) / (2 * self.epsilon)  # V·m²

        # Convert to eV/(10¹³ cm⁻²) using same logic as M1
        slope_display = slope_V_m2 * 1e17  # eV/(10¹³ cm⁻²)

        return slope_display

    def __str__(self):
        return f"M3-Parabolic(εᵣ={self.epsilon_r}, W={self.W_nm:.1f}nm)"
