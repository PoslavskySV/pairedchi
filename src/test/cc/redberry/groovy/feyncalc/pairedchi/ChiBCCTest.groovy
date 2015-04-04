/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2015:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.groovy.Redberry
import org.junit.Test

import static cc.redberry.core.utils.TensorUtils.getAllDiffSimpleTensors
import static cc.redberry.groovy.RedberryStatic.EliminateMetrics
import static cc.redberry.groovy.RedberryStatic.ExpandAll

/**
 * Created by poslavsky on 03/04/15.
 */
class ChiBCCTest {

    @Test
    public void testBottomLeg() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();

            def Mc = '0'.t
            def gcc = 'Vcc_iI = cu[p1_m[charm]]*G_i*T_I*v[p2_m[charm]]'.t
            // (1,2,3)
            Mc += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_bB*eps2^b[h2]*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,2,1)
            Mc += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_bB*eps2^b[h2]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (1,3,2)
            Mc += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_cC*Vcc^cC*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,1,2)
            Mc += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_aA*eps1^a[h1]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (2,3,1)
            Mc += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_cC*Vcc^cC*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t
            // (2,1,3)
            Mc += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_aA*eps1^a[h1]*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t

            Mc <<= 'pCharm_m = p1_m[charm] + p2_m[charm]'.t

            def masses = 'p2_{f}[bottom]*p2^{f}[bottom] = m[bottom]**2'.t & 'p1_{d}[bottom]*p1^{d}[bottom] = m[bottom]**2'.t
            println 'xxxxxx'
            Mc <<= stp.FeynmanRules & stp.spinSingletProjector['bottom']
            Mc <<= stp.dTraceSimplify & stp.uTrace
            Mc <<= masses


            println Mc.size()
            println Mc[0]
            println Mc[1]
            println Mc[2]
            println getAllDiffSimpleTensors(Mc)

            Mc <<= stp.momentums['bottom']
            // Taylor expansion (scalar meson)
            Mc <<= 'q_i[bottom] = q_i'.t.hold & stp.taylor('q_i') & 'q_i = q_i[bottom]'.t

            println Mc.size()
            Mc <<= stp.mandelstam & ExpandAll[EliminateMetrics & stp.mandelstam] &
                    stp.totalSpinProjector['tensor']

            println Mc.size()
            Mc <<= stp.fullSimplify & stp.massesSubs & stp.mFactor & gcc

            println Mc.size()
            println Mc[0]
            println Mc[1]
        }

    }

    @Test
    public void test1123() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();
            for (def spin in ['scalar', 'axial', 'tensor']) {

                println 'GLUON'
                def t = stp.getGluonDiagrams(spin)
                println t.size()
                println t[0]
                println t[1]
                println t[2]

                println 'QUARK'
                t = stp.getQuarkDiagrams(spin)
                println t.size()
                println t[0]
                println t[1]
                println t[2]
            }
        }
    }

    @Test
    public void test1() throws Exception {
        use(Redberry) {
            Setup stp = new Setup(false, true);
            def tensor = "(-4*mc**2+t2-4*mb**2+u1+u2+s+t1)*p2^{i}[charm]*p2^{f}[charm]*p1_{m}[charm]*p1_{d}[charm]*k1_{g}*k1_{e}*G^{dm'}_{l'}*G^{mb'}_{c'}*G^{go'}_{n'}*G^{ed'}_{e'}*G^{cn'}_{m'}*G_{c}^{c'}_{d'}*G^{ke'}_{f'}*G^{bl'}_{k'}*G_{k}^{p'}_{o'}*G_{b}^{a'}_{b'}*T_{C}^{B'}_{C'}*T^{CK'}_{J'}*T_{A}^{C'}_{D'}*T^{BJ'}_{I'}*T_{B}^{A'}_{B'}*T^{AL'}_{K'}*v^{k'I'}[(1/2)*k2_{h}-(1/2)*p2_{h}[charm]-(1/2)*p1_{h}[charm]+(1/2)*k1_{h}]*v^{f'D'}[(1/2)*k2_{h}-(1/2)*p2_{h}[charm]-(1/2)*p1_{h}[charm]+(1/2)*k1_{h}]*q_{i}[bottom]*q_{f}[bottom]*cu_{a'A'}[(1/2)*k2_{g}-(1/2)*p2_{g}[charm]-(1/2)*p1_{g}[charm]+(1/2)*k1_{g}]*cu_{p'L'}[(1/2)*k2_{g}-(1/2)*p2_{g}[charm]-(1/2)*p1_{g}[charm]+(1/2)*k1_{g}]+2*p2^{h}[charm]*p2^{i}[charm]*p2^{f}[charm]*p1_{d}[charm]*p1_{m}[charm]*p1^{c}[charm]*k1_{e}*k1_{g}*G_{h}^{n'}_{m'}*G_{c}^{c'}_{d'}*G^{dm'}_{l'}*G^{mb'}_{c'}*G^{go'}_{n'}*G^{ed'}_{e'}*G^{bl'}_{k'}*G^{ke'}_{f'}*G_{b}^{a'}_{b'}*G_{k}^{p'}_{o'}*T^{CK'}_{J'}*T_{C}^{B'}_{C'}*T_{A}^{C'}_{D'}*T^{BJ'}_{I'}*T_{B}^{A'}_{B'}*T^{AL'}_{K'}*v^{f'D'}[(1/2)*k2_{h}-(1/2)*p2_{h}[charm]-(1/2)*p1_{h}[charm]+(1/2)*k1_{h}]*v^{k'I'}[(1/2)*k2_{h}-(1/2)*p2_{h}[charm]-(1/2)*p1_{h}[charm]+(1/2)*k1_{h}]*q_{i}[bottom]*q_{f}[bottom]*cu_{a'A'}[(1/2)*k2_{g}-(1/2)*p2_{g}[charm]-(1/2)*p1_{g}[charm]+(1/2)*k1_{g}]*cu_{p'L'}[(1/2)*k2_{g}-(1/2)*p2_{g}[charm]-(1/2)*p1_{g}[charm]+(1/2)*k1_{g}]+2*p2^{c}[charm]*p2^{i}[charm]*p2^{f}[charm]*p1_{d}[charm]*p1_{m}[charm]*p1^{h}[charm]*k1_{e}*k1_{g}*G_{c}^{c'}_{d'}*G_{h}^{n'}_{m'}*G^{dm'}_{l'}*G^{mb'}_{c'}*G^{go'}_{n'}*G^{ed'}_{e'}*G^{bl'}_{k'}*G^{ke'}_{f'}*G_{b}^{a'}_{b'}*G_{k}^{p'}_{o'}*T^{CK'}_{J'}*T_{C}^{B'}_{C'}*T_{A}^{C'}_{D'}*T^{BJ'}_{I'}*T_{B}^{A'}_{B'}*T^{AL'}_{K'}*v^{f'D'}[(1/2)*k2_{h}-(1/2)*p2_{h}[charm]-(1/2)*p1_{h}[charm]+(1/2)*k1_{h}]*v^{k'I'}[(1/2)*k2_{h}-(1/2)*p2_{h}[charm]-(1/2)*p1_{h}[charm]+(1/2)*k1_{h}]*q_{i}[bottom]*q_{f}[bottom]*cu_{a'A'}[(1/2)*k2_{g}-(1/2)*p2_{g}[charm]-(1/2)*p1_{g}[charm]+(1/2)*k1_{g}]*cu_{p'L'}[(1/2)*k2_{g}-(1/2)*p2_{g}[charm]-(1/2)*p1_{g}[charm]+(1/2)*k1_{g}]".t
            tensor <<= stp.epsSum & stp.fullSimplify &
                    stp.uTrace & stp.dTraceSimplify &
                    stp.fullSimplify & stp.massesSubs

            println tensor.indices.size()
            println tensor.class
            println getAllDiffSimpleTensors(tensor)
            println tensor
        }
    }

    @Test
    public void testWardIdentities() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC()
            for (def charmSpin in ['scalar', 'axial', 'tensor'])
                for (def bottomSpin in ['scalar', 'axial', 'tensor']) {
                    def M = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
                    for (def i in [1, 2])
                        assert ({
                            def ward = ("eps${i}_i[h${i}] = k${i}_i".t
                                    & "k${3 - i}_i = p_i[bottom] - k${i}_i + p_i[charm]".t
                            ) >> M
                            ward <<= stp.fullSimplify & stp.massesSubs & stp.mFactor
                            ward
                        }()) == 0.t
                }
        }
    }

    @Test
    public void testCalculateAll() throws Exception {
        def stp = new ChiBCC()
        stp.calculateAll()
    }

    @Test
    public void test2() throws Exception {
        use(Redberry) {
            def stp = new ChiBCC()
            def mc = stp.mc('axial')

            println mc
        }
    }
}
