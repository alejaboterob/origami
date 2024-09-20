# -*- coding: utf-8 -*-
"""
Written by: Yasuhiro Miyazawa
"""

# Import packages
import numpy as np
import numpy.typing as npt
from scipy import optimize
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt

import sys
import os
import shutil
import zipfile
import glob

from tqdm import tqdm

# Figure parameters
# plt.style.use('./common/custom.mplstyle')


class ReschOrigamiAnalysis:
    def __init__(self, dir_save, fname_vtk):
        print('{0:*^49}'.format(''))
        print('* {0:^45} *'.format('Resch Origami Static Analysis Tool'))
        print('{0:*^49}'.format(''))

        # axis 60 rotation about z on hexagon
        theta1 = np.pi / 3
        T_4 = np.array([[np.cos(theta1), -np.sin(theta1), 0], [np.sin(theta1), np.cos(theta1), 0], [0, 0, 1]])
        # axis 120 rotation about z on hexagon
        theta2 = 2 * np.pi / 3
        T_5 = np.array([[np.cos(theta2), np.sin(theta2), 0], [-np.sin(theta2), np.cos(theta2), 0], [0, 0, 1]])
        T_6 = T_4.T

        self.T_4 = T_4
        self.T_5 = T_5
        self.T_6 = T_6

        self.dir_save = dir_save
        self.fname_vtk = fname_vtk

        pass

    def geo_init(self, La: float,
                 theta_M0: float,
                 theta_MKJS: npt.ArrayLike,
                 use_lookuptab: bool = True):
        n_orbit = 1
        n_hex = 1 + 3 * (n_orbit + n_orbit**2)
        n_node_hex = int(6 * n_hex)
        n_node_hex_origin = int(n_hex)
        n_node_tri = int(6 * (n_orbit + n_orbit**2))
        n_node = n_node_hex + n_node_tri  # + n_node_hex_origin
        n_edge_hex = int(6 * n_hex)
        n_edge_tri = int(7 * n_node_tri)
        n_edge = n_edge_hex + n_edge_tri
        n_poly_hex = n_hex
        n_poly_tri = int(6 * (n_node_tri - n_orbit * 3))
        n_poly = n_poly_hex + n_poly_tri

        # Determine crease length
        Lb = La / np.sqrt(3.)
        Lc = 2. * La / np.sqrt(3.)

        self.La = La
        self.Lb = Lb
        self.Lc = Lc
        self.n_orbit = n_orbit
        self.theta_M0 = theta_M0
        #
        self.n_node_hex = n_node_hex
        self.n_node_tri = n_node_tri
        self.n_node = n_node
        #
        self.n_edge_hex = n_edge_hex
        self.n_edge_tri = n_edge_tri
        self.n_edge = n_edge
        #
        self.n_poly_hex = n_poly_hex
        self.n_poly_tri = n_poly_tri
        self.n_poly = n_poly

        # Import fold angles
        if use_lookuptab:
            data_th = np.genfromtxt('./lookuptable/theta_MKJS.csv', delimiter=',')
            theta_M = data_th[:, 0]
            theta_K = data_th[:, 1]
            theta_J = data_th[:, 2]
            theta_S = data_th[:, 3]
        else:
            theta_M = theta_MKJS[:, 0]
            theta_K = theta_MKJS[:, 1]
            theta_J = theta_MKJS[:, 2]
            theta_S = theta_MKJS[:, 3]

        self.f_MK = interpolate.interp1d(theta_M, theta_K, kind='cubic')
        self.f_MJ = interpolate.interp1d(theta_M, theta_J, kind='cubic')
        self.f_MS = interpolate.interp1d(theta_M, theta_S, kind='cubic')

        theta_K0 = self.f_MK(theta_M0)
        theta_J0 = self.f_MJ(theta_M0)
        theta_S0 = self.f_MS(theta_M0)

        if n_orbit == 1:
            vert_xyz, EdgeConct, Polyg_tri, Polyg_hex,\
                PolygAdj, Polyg2Edge = self.geo_init_1orbit(La, theta_M0)

        print('{0:*^49}'.format(''))
        print(' Initialize Resch tessellation')
        print('  {0:<32s} : {1:<12d}'.format('Number of orbit', n_orbit))
        print('  {0:<32s} : {1:<.12f}'.format('Hexagon side length (a)', La))
        print('  {0:<32s} : {1:<.12f}'.format('Minor crease length (b)', Lb))
        print('  {0:<32s} : {1:<.12f}'.format('Major crease length (c)', Lc))
        print('  {0:<32s}'.format('Initial fold angles th_M0'))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_M0', np.degrees(theta_M0)))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_K0', np.degrees(theta_K0)))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_J0', np.degrees(theta_J0)))
        print('    {0:<30s} : {1:<.8f} (deg)'.format('th_S0', np.degrees(theta_S0)))
        print('  {0:<32s}'.format('Node, crease, polygon'))
        print('    {0:<30s} : {1:d}, {2:d}, {3:d}'.format('Total', n_node, n_edge, n_poly))
        print('      {0:<28s} : {1:d}, {2:d}, {3:d}'.format('Hexagon', n_node_hex, n_edge_hex, n_poly_hex))
        print('      {0:<28s} : {1:d}, {2:d}, {3:d}'.format('Triangle', n_node_tri, n_edge_tri, n_poly_tri))
        print('{0:*^49}'.format(''))

        return

    def geo_init_1orbit(self, La, theta_M0):

        self.La = La
        # self.n_node = 54

        # n_tri = 48  # 6*6+2*6=48
        # n_hex = 7
        n_tri = self.n_poly_tri
        n_hex = self.n_poly_hex

        EdgeConct = np.zeros((self.n_edge, 2), dtype=int)
        PolygAdj = np.zeros((self.n_edge, 2), dtype=int)
        Polyg2Edge = np.zeros((self.n_poly, 2), dtype=int)

        # Define triangles
        Polyg_tri = np.zeros((n_tri, 3), dtype=int)

        Polyg_tri[0, :] = [0, 1, 43]
        Polyg_tri[1, :] = [1, 2, 44]
        Polyg_tri[2, :] = [2, 3, 45]
        Polyg_tri[3, :] = [3, 4, 46]
        Polyg_tri[4, :] = [4, 5, 47]
        Polyg_tri[5, :] = [5, 0, 42]
        Polyg_tri[6, :] = [0, 9, 42]
        Polyg_tri[7, :] = [1, 16, 43]
        Polyg_tri[8, :] = [2, 23, 44]
        Polyg_tri[9, :] = [3, 24, 45]
        Polyg_tri[10, :] = [4, 31, 46]
        Polyg_tri[11, :] = [5, 38, 47]
        Polyg_tri[12, :] = [0, 9, 43]
        Polyg_tri[13, :] = [1, 16, 44]
        Polyg_tri[14, :] = [2, 23, 45]
        Polyg_tri[15, :] = [3, 24, 46]
        Polyg_tri[16, :] = [4, 31, 47]
        Polyg_tri[17, :] = [5, 38, 42]
        Polyg_tri[18, :] = [8, 9, 43]
        Polyg_tri[19, :] = [15, 16, 44]
        Polyg_tri[20, :] = [22, 23, 45]
        Polyg_tri[21, :] = [24, 29, 46]
        Polyg_tri[22, :] = [30, 31, 47]
        Polyg_tri[23, :] = [37, 38, 42]
        Polyg_tri[24, :] = [8, 17, 43]
        Polyg_tri[25, :] = [16, 17, 43]
        Polyg_tri[26, :] = [15, 18, 44]
        Polyg_tri[27, :] = [18, 23, 44]
        Polyg_tri[28, :] = [22, 25, 45]
        Polyg_tri[29, :] = [24, 25, 45]
        Polyg_tri[30, :] = [29, 32, 46]
        Polyg_tri[31, :] = [31, 32, 46]
        Polyg_tri[32, :] = [30, 39, 47]
        Polyg_tri[33, :] = [38, 39, 47]
        Polyg_tri[34, :] = [10, 37, 42]
        Polyg_tri[35, :] = [9, 10, 42]
        Polyg_tri[36, :] = [8, 17, 48]
        Polyg_tri[37, :] = [12, 17, 48]
        Polyg_tri[38, :] = [15, 18, 49]
        Polyg_tri[39, :] = [18, 19, 49]
        Polyg_tri[40, :] = [25, 26, 50]
        Polyg_tri[41, :] = [22, 25, 50]
        Polyg_tri[42, :] = [32, 33, 51]
        Polyg_tri[43, :] = [29, 32, 51]
        Polyg_tri[44, :] = [30, 39, 52]
        Polyg_tri[45, :] = [39, 40, 52]
        Polyg_tri[46, :] = [10, 11, 53]
        Polyg_tri[47, :] = [10, 37, 53]

        # Define hexagons
        Polyg_hex = np.zeros((n_hex, 6), dtype=int)
        for i in range(self.n_poly_hex):
            Polyg_hex[i, :] = np.arange(6 * i, 6 * (i + 1), dtype=int)

        # Define edge connectivity
        iter = 0
        # Hexagon
        for ip in range(self.n_poly_hex):
            for ii in range(6):
                EdgeConct[iter, :] = np.sort([Polyg_hex[ip, ii - 1], Polyg_hex[ip, ii]])
                iter = iter + 1

        # Triangle
        for ip in range(self.n_poly_tri):
            for ii in range(3):
                idx = np.sort([Polyg_tri[ip, ii - 1], Polyg_tri[ip, ii]])
                idx_test = np.prod(np.isin(EdgeConct[:iter, :], idx), axis=1, dtype=bool)
                if idx_test.any():
                    pass
                else:
                    EdgeConct[iter, :] = idx
                    iter = iter + 1

        self.Polyg_tri = Polyg_tri
        self.Polyg_hex = Polyg_hex
        self.n_poly_tri = n_tri
        self.n_poly_hex = n_hex

        self.EdgeConct = EdgeConct

        theta_M = np.array([theta_M0, ])
        vert_xyz = self.solve_geo_1orbit(theta_M=theta_M, vtk_out=False, fig_out=False)

        self.Polyg = np.zeros((self.n_poly), dtype=object)
        for i in range(self.n_poly_hex):
            self.Polyg[i] = Polyg_hex[i, :]
        for i in range(self.n_poly_tri):
            self.Polyg[self.n_poly_hex + i] = Polyg_tri[i, :]

        return vert_xyz, EdgeConct, Polyg_tri, Polyg_hex, PolygAdj, Polyg2Edge

    def calc_ukvector(self, theta_K: float, theta_M: float, La: float):
        u1_b = La * np.array([0.5 * np.sin(theta_M), 0.5 / np.sqrt(3.), 0.5 * np.cos(theta_M)])
        u2_b = La * np.array([-0.5 * np.sin(theta_M), 0.5 * np.sqrt(3.), -0.5 * np.cos(theta_M)])

        n1 = u2_b / La

        mu1 = 2. * theta_K
        k1_b = np.dot(u1_b, n1) * n1 * (1. - np.cos(mu1)) + np.cos(mu1) * u1_b - np.sin(mu1) * np.cross(n1, u1_b)

        E = u1_b + u2_b - k1_b
        G = u1_b
        EG = G - E
        # n2 = EG / (2. * a / np.sqrt(3.))
        n2 = 0.5 * EG * np.sqrt(3) / La
        mu2 = 2 * theta_M
        k4_b = np.dot(k1_b, n2) * n2 * (1 - np.cos(mu2)) + np.cos(mu2) * k1_b - np.sin(mu2) * np.cross(n2, k1_b)
        k3_b = -(u2_b - k1_b + k4_b)
        u3_b = La * np.array([0.5 * np.sin(theta_M), 0.5 * np.sqrt(3.), -0.5 * np.cos(theta_M)])

        psi = np.pi - np.arccos(np.sin(theta_M)) - np.arccos(1. / (6. * np.sqrt(1. / 9. - 1. / 12. * np.sin(theta_M)**2))) - \
            np.emath.arccos((2. / 3. - np.sin(theta_M)**2) / (2. * np.sqrt(1. / 9. - 1. / 12. * np.sin(theta_M)**2))).real

        T_1 = np.array([[1., 0., 0.], [0., np.cos(psi), -np.sin(psi)], [0., np.sin(psi), np.cos(psi)]])
        t_1 = T_1.T
        p1 = 2. * np.pi / 3.
        T_3 = np.array([[np.cos(p1), np.sin(p1), 0.], [-np.sin(p1), np.cos(p1), 0.], [0., 0., 1.]])
        u4_b = np.dot(np.dot(np.dot(t_1, T_3), T_1), u3_b)

        return psi, u1_b, u2_b, u3_b, u4_b, k1_b, k3_b, k4_b

    def funcMK(self, theta_K: float, theta_M: float, La: float):

        _, _, _, _, u4_b, _, k3_b, _ = self.calc_ukvector(theta_K, theta_M, La)

        f = np.dot(u4_b, k3_b) - La**2 * np.cos(np.pi / 3.)

        return f

    def calc_theta_JS(self, theta_M: float, theta_K: float, La: float):

        psi, u1_b, _, _, u4_b, _, k3_b, _ = self.calc_ukvector(theta_K, theta_M, La)

        v1_b = k3_b - 0.5 * u4_b
        theta_J = np.emath.arccos(2. * np.dot(u1_b, v1_b) / La**2).real
        theta_S = np.arcsin(np.cos(psi))

        return theta_J, theta_S

    def find_theta_MKJS(self, theta_M: list | np.ndarray,
                        La: float = 1.0,
                        x0: float = 0.0,
                        method: str = 'hybr',
                        fig_out: bool = True,
                        save_lookuptab: bool = False
                        ):
        '''
        find_theta_MKJS
        Determines theta_K and theta_J as a function of theta_M

        Args:
            theta_M (list | np.ndarray): Array of theta_M
            a (float, optional): Side length. Defaults to 1.0.
            x0 (float, optional): Initial guess for the first step. Defaults to 0.0.
            method (str, optional): Method for scipy.optimize.root. Defaults to 'hybr'.
            fig_out (bool, optional): If True, displays figures. Defaults to True.
            save_lookuptab (bool, optional): If True, saves lookup table. Defaults to False.

        Returns:
            (list | np.ndarray): List of theta_K and theta_J
        '''

        niter = len(theta_M)

        # Allocate arrays to store solutions
        theta_J = np.zeros_like(theta_M)
        theta_K = np.zeros_like(theta_M)
        theta_S = np.zeros_like(theta_M)

        for i in tqdm(range(niter)):
            # Find root of funcMK at each theta_M
            output = optimize.root(self.funcMK, x0, args=(theta_M[i], La),
                                   method=method)
            # Save solution
            theta_K[i] = output.x
            # Update initial guess for next step
            x0 = theta_K[i]

        # Substitute first and last entry with correct mathematical values
        theta_K[0] = 0.0
        theta_K[-1] = 0.5 * np.pi

        # Calculate theta_J and theta_S
        for i in range(niter):
            theta_J[i], theta_S[i] = self.calc_theta_JS(theta_M=theta_M[i], theta_K=theta_K[i], La=La)

        # Store fold angles
        theta_MKJS = np.zeros((4, len(theta_M)))
        theta_MKJS[0, :] = theta_M
        theta_MKJS[1, :] = theta_K
        theta_MKJS[2, :] = theta_J
        theta_MKJS[3, :] = theta_S

        # Save lookup table
        if save_lookuptab:
            print(' Save lookup table')
            fpath = './lookuptable/theta_MKJS.csv'
            file_exists = os.path.exists(fpath)
            if file_exists:
                with open(fpath, 'r') as fp:
                    for count, line in enumerate(fp):
                        pass
                count = count + 1
                file_stat = os.stat(fpath)
                if file_stat.st_size >= 1e3 and file_stat.st_size < 1e6:
                    fsize = file_stat.st_size * 1e-3
                    funit = 'KB'
                elif file_stat.st_size >= 1e6:
                    fsize = file_stat.st_size * 1e-6
                    funit = 'MB'
                else:
                    fsize = file_stat.st_size
                    funit = 'B'
                print('  ! File already exists ({0:d} lines; {1:.2f} {2:2s}) !'.format(count, fsize, funit))
                value = input('  Overwrite ? [Press Y to overwrite, press elsewhere otherwise]\n')
                if value == 'Y' or value == 'y':
                    print('   Saving lookup table...')
                    data = np.stack((theta_M, theta_K, theta_J, theta_S), axis=-1)
                    np.savetxt(fpath, data, delimiter=',')
                    print('    ...Saved')
                else:
                    print('   File not saved.')
                    pass
            else:
                print('   Saving lookup table...')
                data = np.stack((theta_M, theta_K, theta_J, theta_S), axis=-1)
                np.savetxt(fpath, data, delimiter=',')
                print('    ...Saved')

        # Generate figures
        if fig_out:
            plt.figure('theta MKJS')
            plt.plot(np.degrees(theta_M), np.degrees(theta_K), label='$\\theta_K$')
            plt.plot(np.degrees(theta_M), np.degrees(theta_J), label='$\\theta_J$')
            plt.plot(np.degrees(theta_M), np.degrees(theta_S), label='$\\theta_S$')
            plt.axhline(y=90.0, linewidth=1.0, color='#484848', linestyle='dashed')
            plt.legend(loc=4)
            plt.xlim(np.degrees(min(theta_M)), np.degrees(max(theta_M)))
            plt.ylim(0, 180)
            plt.xlabel('$\\theta_M$ (deg)')
            plt.ylabel('$\\theta_K$, $\\theta_J$ (deg)')

        return theta_MKJS

    def BaseHexLayer(self, theta_M, theta_J):
        A = np.zeros((8, 6, 3))
        B = np.zeros((8, 6, 3))
        C = np.zeros((8, 6, 3))
        D = np.zeros((8, 6, 3))
        E = np.zeros((8, 6, 3))

        A[0, 1, :] = np.array([0, self.La, 0])
        A[0, 0, :] = np.array([-0.5 * np.sqrt(3) * self.La, 0.5 * self.La, 0])
        A[0, 5, :] = np.array([-0.5 * np.sqrt(3) * self.La, -0.5 * self.La, 0])

        A1B1 = self.La / np.sqrt(3) * np.array([1 / 2 * np.cos(theta_J), -np.sqrt(3) / 2 * np.cos(theta_J), -np.sin(theta_J)])
        B[0, 1, :] = A[0, 1, :] + A1B1
        B1A1 = -A1B1
        A0B1 = B[0, 1] - A[0, 0]
        n3 = A0B1 / np.linalg.norm(A0B1)
        mu3 = 2 * theta_M
        B1C1 = np.dot(B1A1, n3) * n3 * (1 - np.cos(mu3)) + np.cos(mu3) * B1A1 - np.sin(mu3) * np.cross(n3, B1A1)
        C[0, 1] = B[0, 1] + B1C1
        C1A0 = A[0, 0] - C[0, 1]
        n4 = B1C1 / np.linalg.norm(B1C1)
        psi = np.pi - np.arccos(np.sin(theta_M)) - np.arccos(1. / (6 * (1 / 9 - 1 / 12 * np.sin(theta_M)**2)**0.5)) - \
            np.emath.arccos((2 / 3 - np.sin(theta_M)**2) / (2 * (1 / 9 - 1 / 12 * np.sin(theta_M)**2)**0.5)).real
        thetaS = np.arcsin(np.cos(psi))
        mu4 = 2 * thetaS
        C1E1 = np.dot(C1A0, n4) * n4 * (1 - np.cos(mu4)) + np.cos(mu4) * C1A0 - np.sin(mu4) * np.cross(n4, C1A0)
        E[0, 1] = C[0, 1] + C1E1
        E1B1 = B[0, 1] - E[0, 1]
        n5 = E1B1 / np.linalg.norm(E1B1)
        mu5 = mu3
        B1D1 = np.dot(B1C1, n5) * n5 * (1 - np.cos(mu5)) + np.cos(mu5) * B1C1 - np.sin(mu5) * np.cross(n5, B1C1)
        D[0, 1] = B[0, 1] + B1D1

        B[0, 0] = self.T_4 @ B[0, 1]
        C[0, 0] = self.T_4 @ C[0, 1]
        E[0, 0] = self.T_4 @ E[0, 1]
        D[0, 0] = self.T_4 @ D[0, 1]
        A[0, 2] = self.T_5 @ A[0, 0]
        B[0, 2] = self.T_5 @ B[0, 0]
        C[0, 2] = self.T_5 @ C[0, 0]
        E[0, 2] = self.T_5 @ E[0, 0]
        D[0, 2] = self.T_5 @ D[0, 0]
        A[0, 3] = self.T_5 @ A[0, 1]
        B[0, 3] = self.T_5 @ B[0, 1]
        C[0, 3] = self.T_5 @ C[0, 1]
        E[0, 3] = self.T_5 @ E[0, 1]
        D[0, 3] = self.T_5 @ D[0, 1]
        A[0, 4] = self.T_5 @ A[0, 2]
        B[0, 4] = self.T_5 @ B[0, 2]
        C[0, 4] = self.T_5 @ C[0, 2]
        E[0, 4] = self.T_5 @ E[0, 2]
        D[0, 4] = self.T_5 @ D[0, 2]
        A[0, 5] = self.T_5 @ A[0, 3]
        B[0, 5] = self.T_5 @ B[0, 3]
        C[0, 5] = self.T_5 @ C[0, 3]
        E[0, 5] = self.T_5 @ E[0, 3]
        D[0, 5] = self.T_5 @ D[0, 3]

        return A, B, C, D, E

    def M01_1(self, angles, A, C, D, E):
        r1 = np.array([[1, 0, 0], [0, np.cos(angles[0]), np.sin(angles[0])], [0, -np.sin(angles[0]), np.cos(angles[0])]])
        r2 = np.array([[np.cos(angles[1]), 0, -np.sin(angles[1])], [0, 1, 0], [np.sin(angles[1]), 0, np.cos(angles[1])]])
        r3 = r2 @ r1

        N1 = np.cross((A[0, 2] - A[0, 3]), (A[0, 4] - A[0, 3]))
        N2 = np.cross((E[0, 1] - C[0, 1]), (D[0, 0] - C[0, 1]))
        f = r3 @ N1
        f = f - N2

        return f[0:2]

    def M01_2(self, angles, Rangles, A, C, E):
        R1 = np.array([[1, 0, 0], [0, np.cos(Rangles[0]), np.sin(Rangles[0])], [0, -np.sin(Rangles[0]), np.cos(Rangles[0])]])
        R2 = np.array([[np.cos(Rangles[1]), 0, -np.sin(Rangles[1])], [0, 1, 0], [np.sin(Rangles[1]), 0, np.cos(Rangles[1])]])
        r3 = np.array([[np.cos(angles[0]), -np.sin(angles[0]), 0], [np.sin(angles[0]), np.cos(angles[0]), 0], [0, 0, 1]])
        r4 = r3 @ R2 @ R1

        N1 = (A[0, 2] - A[0, 3])
        N2 = (E[0, 1] - C[0, 1])
        f = r4 @ N1
        f = f - N2

        return f[1]

    def M01(self, x0, theta_M, theta_J):
        A, B, C, D, E = self.BaseHexLayer(theta_M, theta_J)

        Rangles = np.zeros(3)

        sol = optimize.root(self.M01_1, x0[0:2], args=(A, C, D, E))
        Rangles[0] = sol.x[0]
        Rangles[1] = sol.x[1]

        sol = optimize.root(self.M01_2, x0[2], args=(Rangles, A, C, E))
        Rangles[2] = sol.x

        R1 = np.array([[1, 0, 0], [0, np.cos(Rangles[0]), np.sin(Rangles[0])], [0, -np.sin(Rangles[0]), np.cos(Rangles[0])]])  # ccw
        R2 = np.array([[np.cos(Rangles[1]), 0, -np.sin(Rangles[1])], [0, 1, 0], [np.sin(Rangles[1]), 0, np.cos(Rangles[1])]])  # ccw
        R3 = np.array([[np.cos(Rangles[2]), -np.sin(Rangles[2]), 0], [np.sin(Rangles[2]), np.cos(Rangles[2]), 0], [0, 0, 1]])  # cw

        matR = R3 @ R2 @ R1

        mat_pack = [A, B, C, D, E]

        return matR, Rangles, mat_pack  # 0th-1st hexagon radial transformation

    def solve_geo_1orbit(self, theta_M: npt.ArrayLike,
                         vtk_out: bool = False, fig_out: bool = True,
                         progbar_out: bool = True,
                         save_zip: bool = True):

        progbar_disable = ~progbar_out

        niter = len(theta_M)

        theta_J = self.f_MJ(theta_M)

        A_global = np.zeros((42, 3))
        B_global = np.zeros((42, 3))

        vert_xyz = np.zeros((niter, self.n_node, 3))

        x0 = [0., 0., 0.]

        for i in tqdm(range(niter), disable=progbar_disable):
            matR, Rangles, mat_pack = self.M01(x0, theta_M[i], theta_J[i])
            A, B, C, D, E = mat_pack

            # Update initial guess
            # x0 = Rangles

            for ii in range(0, 6):
                A[1, ii] = E[0, 1] + matR @ (A[0, ii] - A[0, 2])
                B[1, ii] = E[0, 1] + matR @ (B[0, ii] - A[0, 2])
                C[1, ii] = E[0, 1] + matR @ (C[0, ii] - A[0, 2])
                D[1, ii] = E[0, 1] + matR @ (D[0, ii] - A[0, 2])
                E[1, ii] = E[0, 1] + matR @ (E[0, ii] - A[0, 2])
            for ii in range(1, 6):
                for jj in range(0, 5):
                    A[ii + 1, jj + 1] = self.T_6 @ A[ii, jj]
                    A[ii + 1, 0] = self.T_6 @ A[ii, 5]
                    B[ii + 1, jj + 1] = self.T_6 @ B[ii, jj]
                    B[ii + 1, 0] = self.T_6 @ B[ii, 5]
                    C[ii + 1, jj + 1] = self.T_6 @ C[ii, jj]
                    C[ii + 1, 0] = self.T_6 @ C[ii, 5]
                    D[ii + 1, jj + 1] = self.T_6 @ D[ii, jj]
                    D[ii + 1, 0] = self.T_6 @ D[ii, 5]
                    E[ii + 1, jj + 1] = self.T_6 @ E[ii, jj]
                    E[ii + 1, 0] = self.T_6 @ E[ii, 5]

            for ii in range(0, 7):
                for jj in range(0, 6):
                    A_global[jj + 6 * ii] = A[ii, jj]
                    B_global[jj + 6 * ii] = B[ii, jj]

            vert_xyz[i, 0:42, :] = A_global
            vert_xyz[i, 42:48, :] = B_global[0:6]
            vert_xyz[i, 48:54, :] = B_global[[8, 15, 22, 29, 30, 37]]

        if vtk_out:
            if os.path.exists(self.dir_save):
                shutil.rmtree(self.dir_save)
            os.makedirs(self.dir_save)
            print('Create VTK...')
            for i in tqdm(range(niter)):
                self.write_vtk_1orbit(i, vert_xyz[i, :, :])

            if save_zip:
                zp = zipfile.ZipFile('%s.zip' % self.dir_save, 'w')
                dfile = glob.glob('%s/*.vtk' % self.dir_save)
                dfile = np.sort(dfile)
                for i in range(len(dfile)):
                    zp.write(filename=dfile[i], arcname=None, compress_type=None, compresslevel=9)
                zp.close()

        if fig_out:
            from plot_geometry import plot_projection
            plot_projection(vert_xyz=vert_xyz[0],
                            Polyg=self.Polyg,
                            EdgeConct=self.EdgeConct,
                            figname='(Resck63O1-init)')
            plot_projection(vert_xyz=vert_xyz[-1],
                            Polyg=self.Polyg,
                            EdgeConct=self.EdgeConct,
                            figname='(Resck63O1-final)')

        return vert_xyz

    def write_vtk_1orbit(self, fnum: int, vert_xyz: list):  # , strain: list ):

        with open('%s/%s_%05d.vtk' % (self.dir_save, self.fname_vtk, fnum), 'w') as f:
            f.write('# vtk DataFile Version 3.0\n')
            f.write('Hexagon-triangle Resch origami\n')
            f.write('ASCII\n')
            num_points = np.size(vert_xyz, axis=0)
            f.write('DATASET POLYDATA\n')
            f.write('POINTS %d double\n' % num_points)
            for i in range(num_points):
                f.write("%.15e %.15e %.15e\n" % (vert_xyz[i, 0], vert_xyz[i, 1], vert_xyz[i, 2]))

            num_dataset = self.n_poly_tri + self.n_poly_hex
            num_datanum = 4 * self.n_poly_tri + 7 * self.n_poly_hex
            f.write('POLYGONS %d %d\n' % (num_dataset, num_datanum))
            # POLYGONS:triangles
            for ip in range(self.n_poly_tri):
                f.write('3 %d %d %d\n' % (self.Polyg_tri[ip, 0], self.Polyg_tri[ip, 1], self.Polyg_tri[ip, 2]))

            # POLYGONS:hexagons
            for ip in range(self.n_poly_hex):
                f.write('6 %d %d %d %d %d %d\n' % (self.Polyg_hex[ip, 0], self.Polyg_hex[ip, 1],
                        self.Polyg_hex[ip, 2], self.Polyg_hex[ip, 3], self.Polyg_hex[ip, 4], self.Polyg_hex[ip, 5]))

        return


if __name__ == "__main__":

    theta0 = np.radians(0.0)
    thetaf = np.radians(90.0)
    nbin = 64

    theta_M = np.linspace(theta0, thetaf, nbin)

    h0 = 1.0
    ph0 = 60.0
    n_unit = 1
    n_side = 6
    La = 1.0    

    dir_save = './vtk_kresling'
    fname_vtk = 'kresling'



    roa = ReschOrigamiAnalysis(dir_save, fname_vtk)

    xx = np.linspace(0.0, 0.5 * np.pi, int(2 * nbin))
    theta_MKJS = roa.find_theta_MKJS(theta_M=xx,
                                            La=La, x0=np.radians(0.1),
                                            save_lookuptab=True,
                                            )
    roa.geo_init(La=La,
                theta_M0=theta_M[0],
                theta_MKJS=theta_MKJS.T,
                use_lookuptab=False,
                )

    roa.solve_geo_1orbit(theta_M=theta_M,
                            vtk_out=False,
                            fig_out=True,
                            save_zip=False,
                            )




# import json
# import numpy as np

# from plot_ori_panels import plot_panels, save_gif_PIL
# import matplotlib.pyplot as plt

# from ConfigYoshimura import plot_yoshimura_nodes

# def ConfigResch(file_path):
#     try:
#         # Open and read the file as a regular text file
#         with open(file_path, 'r') as file:
#             file_content = file.read()  # Read the entire content of the file

#         # Attempt to parse the content as JSON
#         try:
#             data_dict = json.loads(file_content)  # This loads the JSON content into a dictionary
#             print(data_dict)  # Output the dictionary to confirm it's loaded correctly
#         except json.JSONDecodeError:
#             print("Error decoding JSON. The file content might not be in JSON format.")

#     except FileNotFoundError:
#         print(f"File not found: {file_path}")

#     Nodes = np.array(data_dict['vertices_coords'])
#     Panels = np.array(data_dict['faces_vertices'])
#     fda = np.array(data_dict['edges_foldAngle']) 
#     BDRY = np.where(np.array(fda) == None)

#     return Nodes, Panels, BDRY




# file_path = r"python\fold\reschTriTessellation _ 20PercentFolded.fold"

# Nodes, Panels, BDRY = ConfigResch(file_path)
# print(Nodes)
# print(Panels)


# try:
#     # Open and read the file as a regular text file
#     with open(file_path, 'r') as file:
#         file_content = file.read()  # Read the entire content of the file

#     # Print the content of the file
#     print(file_content)  # Output the content to confirm it's loaded correctly

# except FileNotFoundError:
#     print(f"File not found: {file_path}")

# plot_yoshimura_nodes(Nodes)

# plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
#                 "facecolors": "#d86a96"}

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# plot_panels(Nodes, Panels, ax, **plot_kwargs)

# mins = np.minimum(Nodes.min(axis=0), Nodes.min(axis=0))
# maxs = np.maximum(Nodes.max(axis=0), Nodes.max(axis=0))

# # plt.xticks([])
# # plt.yticks([])
# ax.auto_scale_xyz([mins[0], maxs[0]],
#             [mins[1], maxs[1]],
#             [mins[2], maxs[2]])    
# # ax.axis('equal')
# # ax.set(xticklabels=[], yticklabels=[], zticklabels=[])


# plt.show(block = False)

# plt.show()

# # import json

# # # Correct file path where your file is located
# # file_path = 'C:/Users/mboter45/OneDrive - Universidad EAFIT/Estructuras origami/fold/waterbomb_10PercentFolded.fold'

# # try:
# #     # Open and read the file as a regular text file
# #     with open(file_path, 'r') as file:
# #         file_content = file.read()  # Read the entire content of the file

# #     # Attempt to parse the content as JSON
# #     try:
# #         data_dict = json.loads(file_content)  # This loads the JSON content into a dictionary
# #         print(data_dict)  # Output the dictionary to confirm it's loaded correctly
# #     except json.JSONDecodeError:
# #         print("Error decoding JSON. The file content might not be in JSON format.")

# # except FileNotFoundError:
# #     print(f"File not found: {file_path}")
