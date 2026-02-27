"""
M2: Fang-Howard variational model (self-consistent field reduction).

Physical picture: Electric field decays with depth, electron cloud width adapts.
Wavefunction: ψ(z) = √(b³/2) · z · exp(-bz/2)
"""

import numpy as np
from physics.constants import Q, EPS0, HBAR
from physics.units import nm_to_m, eV_to_J, J_to_eV, m_to_nm


class FangHowardModel:
    """M2: Fang-Howard variational model with self-consistent iteration."""

    def __init__(self, m_star, epsilon_r, W_nm_initial):
        """
        Initialize the Fang-Howard model.

        Parameters
        ----------
        m_star : float
            Effective mass in kg
        epsilon_r : float
            Relative permittivity
        W_nm_initial : float
            Initial guess for depletion width in nanometers
        """
        self.m_star = m_star
        self.epsilon_r = epsilon_r
        self.W_initial = nm_to_m(W_nm_initial)
        self.epsilon = EPS0 * epsilon_r

        # Store last calculated effective width
        self.W_eff = self.W_initial

    def calculate_b_parameter(self, Es):
        """
        Calculate the variational parameter b.

        Parameters
        ----------
        Es : float
            Surface electric field in V/m

        Returns
        -------
        b : float
            Variational parameter in m⁻¹
        """
        # b = (12·m*·q·Es/ℏ²)^(1/3)
        b = (12 * self.m_star * Q * Es / HBAR**2)**(1/3)
        return b

    def solve_self_consistent(self, Phi_s_eV, max_iter=10, tolerance=1e-4):
        """
        Solve self-consistent equations for given surface potential.

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV
        max_iter : int
            Maximum number of iterations
        tolerance : float
            Convergence tolerance for relative change in Es

        Returns
        -------
        Es : float
            Self-consistent surface field in V/m
        W_eff : float
            Effective width in meters
        n_iter : int
            Number of iterations used
        """
        Phi_s_V = Phi_s_eV  # Potential in volts

        # Initial guess: Es = Phi_s / W_initial (like M1)
        Es = Phi_s_V / self.W_initial

        for iteration in range(max_iter):
            # Calculate b parameter
            b = self.calculate_b_parameter(Es)

            # Calculate effective width: W_eff = 6/b
            W_eff = 6.0 / b

            # Update Es: Phi_s ≈ Es * W_eff / 2  =>  Es = 2*Phi_s / W_eff
            Es_new = 2 * Phi_s_V / W_eff

            # Check convergence
            if abs(Es_new - Es) / Es < tolerance:
                self.W_eff = W_eff  # Store for later use
                return Es_new, W_eff, iteration + 1

            Es = Es_new

        # If didn't converge, still return the last values
        self.W_eff = W_eff
        return Es, W_eff, max_iter

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
        # Handle array input
        if isinstance(Phi_s_eV, np.ndarray):
            return np.array([self.calculate_ns(phi) for phi in Phi_s_eV])

        # Solve self-consistent equations
        Es, W_eff, _ = self.solve_self_consistent(Phi_s_eV)

        # Gauss's law: ns = (epsilon * Es) / q
        ns = (self.epsilon * Es) / Q

        return ns

    def calculate_Phi_s(self, ns, initial_guess=0.3):
        """
        Calculate surface potential from sheet density (inverse calculation).

        This requires numerical root finding since the relationship is implicit.

        Parameters
        ----------
        ns : float
            Sheet density in m⁻²
        initial_guess : float
            Initial guess for Phi_s in eV

        Returns
        -------
        Phi_s_eV : float
            Surface potential in eV
        """
        from scipy.optimize import fsolve

        def residual(Phi_s_eV):
            ns_calculated = self.calculate_ns(Phi_s_eV)
            return ns_calculated - ns

        # Use fsolve to find Phi_s that gives the target ns
        Phi_s_eV = fsolve(residual, initial_guess)[0]
        return Phi_s_eV

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
        # Handle array input
        if isinstance(ns, np.ndarray):
            return np.array([self.calculate_Delta_WF(n) for n in ns])

        # First find Phi_s from ns
        Phi_s_eV = self.calculate_Phi_s(ns)

        # Delta_WF = -Phi_s
        return -Phi_s_eV

    def get_potential(self, Phi_s_eV, z_array):
        """
        Get potential profile V(z).

        For Fang-Howard, we use a soft triangular approximation.

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
        # Solve for self-consistent field
        Es, W_eff, _ = self.solve_self_consistent(Phi_s_eV)

        Phi_s_V = Phi_s_eV  # Potential in volts

        # Electric potential: V_elec(z) = -Phi_s * (1 - z/W_eff) for z ≤ W_eff
        # (Negative because we're measuring from bulk, surface is lower)
        z_ratio = z_array / W_eff
        V_elec = -Phi_s_V * (1 - z_ratio)

        # Potential energy for electron: U = -e * V_elec = +e * Phi_s * (1-z/W_eff)
        # (Positive because electron is negatively charged)
        V_z = -Q * V_elec  # This gives positive energy (barrier for electrons)

        # Set potential to 0 beyond W_eff
        V_z[z_array > W_eff] = 0

        return V_z

    def get_electron_density(self, Phi_s_eV, z_array):
        """
        Get electron density distribution n(z).

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV
        z_array : array
            Array of z positions in meters

        Returns
        -------
        n_z : array
            Electron density in m⁻³ at each z position
        """
        # Solve for self-consistent parameters
        Es, W_eff, _ = self.solve_self_consistent(Phi_s_eV)
        b = self.calculate_b_parameter(Es)

        # Wavefunction: ψ(z) = √(b³/2) · z · exp(-bz/2)
        psi_z = np.sqrt(b**3 / 2) * z_array * np.exp(-b * z_array / 2)

        # Electron density: n(z) = ns * |ψ(z)|²
        ns = self.calculate_ns(Phi_s_eV)
        n_z = ns * psi_z**2

        return n_z

    def get_W_eff_nm(self, Phi_s_eV):
        """
        Get the effective width for a given surface potential.

        Parameters
        ----------
        Phi_s_eV : float
            Surface potential in eV

        Returns
        -------
        W_eff_nm : float
            Effective width in nanometers
        """
        _, W_eff, _ = self.solve_self_consistent(Phi_s_eV)
        return m_to_nm(W_eff)

    def __str__(self):
        return f"M2-Fang-Howard(εᵣ={self.epsilon_r}, W_eff≈{m_to_nm(self.W_eff):.2f}nm)"
