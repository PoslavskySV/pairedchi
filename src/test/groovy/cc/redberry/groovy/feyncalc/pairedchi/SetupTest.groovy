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

import cc.redberry.core.context.CC
import cc.redberry.core.tensor.FastTensors
import cc.redberry.core.tensor.Product
import cc.redberry.core.tensor.Sum
import cc.redberry.core.tensor.SumBuilder
import cc.redberry.core.transformations.Transformation
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.junit.Before
import org.junit.Test

import static cc.redberry.core.indices.IndexType.Matrix1
import static cc.redberry.core.indices.IndexType.Matrix2
import static cc.redberry.groovy.RedberryPhysics.DiracTrace
import static cc.redberry.groovy.RedberryPhysics.mandelstam
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 02/04/15.
 */
class SetupTest {
    @Before
    public void before() throws Exception {
        CC.reset()
    }

    @Test
    public void testInitialisation() {
        use(Redberry) {
            for (def cc in [true, false]) {
                CC.reset()//clear all previous definitions

                def defs = new Setup(cc)
                def t = 'Tr[V_aA*V_bB]'.t
                t <<= defs.FeynmanRules & defs.dTrace & defs.uTrace
                assert t == '-2*g**2*g_{ab}*g_{BA}'.t
            }

        }
    }

    @Test
    public void testTaylor() {
        use(Redberry) {
            def stp = new Setup(false)
            def t = stp.taylor('q_i') >> 'q_i*k^i'.t
            assert t == 'q_i*k^i'.t
        }
    }

    @Test
    public void testEffectiveVertex1() throws Exception {
        use(Redberry) {
            def stp = new Setup(true)
            def vertices = stp.effectiveQuarkoniaVertices()
            for (fl in ['charm', 'bottom'])
                for (spin in ['scalar', 'axial', 'tensor']) {
                    def t = vertices[spin] >> "A${spin}_{aA bB}[${fl}, k1_i, k2_i]".t
                    def expectedContent = [
                            'eps2^{b}[h2]', 'eps1^{a}[h1]',
                            "eps^{f}[h[$fl]]", "eps^{fg}[h[$fl]]", "m[$fl]", "p_{b}[$fl]",
                            'k1^{i}', 'k2^{i}',
                            'e_{abcd}', 'g_{BA}', 'g', 'g_{ab}']
                            .collect({ x -> x.t.hashCode() }).toSet()
                    def actualContent = TensorUtils.getAllDiffSimpleTensors(t).collect({ x -> x.t.hashCode() }).toSet()
                    assert expectedContent.containsAll(actualContent)
                }
        }
    }

    @Test
    public void testEffectiveVertex2() throws Exception {
        use(Redberry) {
            def stp = new Setup(false)
            def vertex = stp.effectivePairVertex()
            for (fl in ['fl', 'charm', 'bottom']) {
                def t = vertex >> "B_{aA bB}[${fl}, k1_i, k2_i]".t
                def expectedContent = [
                        'm[fl]', 'T_{A}', 'g', 'v[p2_{m}[fl]]',
                        'p1_{c}[fl]', 'cu[p1_{m}[fl]]',
                        'G_{b}', 'k1^{c}', 'k2_{m}', 'g_ab', 'g_AB'].t
                        .collect({ x -> x.t.hashCode() }).toSet()
                def actualContent = TensorUtils.getAllDiffSimpleTensors(t).collect({ x -> x.t.hashCode() }).toSet()
                assert expectedContent.containsAll(actualContent)
            }
        }
    }

    /**
     * Testing Ward identities
     */
    @Test
    public void testEffectiveVerticesWard() {
        use(Redberry) {

            def stp = new Setup(true)
            def vertices = stp.effectiveQuarkoniaVertices()

            for (def spin in ['scalar', 'axial', 'tensor']) {
                def vertex = vertices[spin]

                // Checking that A(g + g -> Q) is symmetric
                for (def f in ['charm', 'bottom'])
                    assert (vertex & ExpandAll) >> "A${spin}_{aA bB}[${f}, k1_m, k2_m] - A${spin}_{bB aA}[${f}, k2_m, k1_m]".t == 0.t

                // A(g + g -> Q) does not contains relative momentum
                def A = vertex >> "A${spin}_{aA bB}[fl, k1_m, k2_m]".t
                A.parentAfterChild { t ->
                    assert !('q_i[fl]'.t % t).exists
                }

                // Ward identities
                for (def temp in [
                        vertex >> "k1^a*A${spin}_{aA bB}[charm, k1_i, k2_i]".t,
                        vertex >> "k2^b*A${spin}_{aA bB}[charm, k1_i, k2_i]".t]) {
                    temp = ('k1_a = p_a[charm] - k2_a'.t
                            & stp.fullSimplify
                            & stp.massesSubs
                            & Factor
                    ) >> temp
                    assert temp == 0.t
                }
            }
        }
    }

    @Test
    public void testEffectivePairVertex1() {
        use(Redberry) {
            def stp = new Setup(false)
            def vertex = stp.effectivePairVertex()

            // Checking that A(g + g -> Q) is symmetric
            for (def f in ['charm', 'bottom'])
                assert (vertex & ExpandAll) >> "B_{aA bB}[${f}, k1_m, k2_m] - B_{bB aA}[${f}, k2_m, k1_m]".t == 0.t

        }
    }

    @Test
    public void testEffectivePairVertexWard() {
        use(Redberry) {
            def stp = new Setup(true)

            //two diagrams
            def g2 = 'cu[p1_m[fl]]*(V_aA*D[p1_m[fl] - k1_m, m[fl]]*V_bB + V_bB*D[p1_m[fl] - k2_m, m[fl]]*V_aA)*v[p2_m[fl]]'.t
            //third diagram (3-gluon)
            def g3 = 'cu[p1_m[fl]]*V^cC*v[p2_m[fl]]*G_c^d[k1_a + k2_a]*V_{aA bB dC}[k1_a, k2_a, -k1_a -k2_a]'.t
            def M = g2 + g3
            // Simplifying
            M <<= stp.FeynmanRules & ExpandAll[EliminateMetrics] & EliminateMetrics &
                    'p1_m[fl]*p1^m[fl] = m[fl]**2'.t & 'p2_m[fl]*p2^m[fl] = m[fl]**2'.t
            def pairVertex = 'B_{aA bB}[fl, k1_m, k2_m]'.t.eq(M)


            def mndlst = setMandelstam([k1_m: '0', k2_m: '0', 'p1_m[charm]': 'm[charm]', 'p2_m[charm]': 'm[charm]'])
            def simplifyMetrics = EliminateMetrics &
                    'eps1^a[h1] * k1_a = 0'.t &
                    'eps2^a[h2] * k2_a = 0'.t &
                    mndlst &
                    'd^i_i = 4'.t & 'd^A_A = 8'.t & "d^i'_i' = 4".t & "d^A'_A' = 3".t
            def fullSimplify = simplifyMetrics &
                    ExpandAll[simplifyMetrics] & simplifyMetrics &
                    stp.leviSimplify &
                    ExpandAll[simplifyMetrics] & simplifyMetrics

            def dTraceSimplify = DiracTrace[[Gamma: 'G_a', Simplifications: fullSimplify]]

            for (def amp in [
                    'k1^a * eps2^b[h2] * B_{aA bB}[charm, k1_m, k2_m]'.t,
                    'eps1^a[h1] * k2^b * B_{aA bB}[charm, k1_m, k2_m]'.t
            ]) {
                amp <<= pairVertex & fullSimplify

                def sum = new SumBuilder()
                for (int i = 0; i < amp.size(); ++i)
                    for (int j = 0; j < amp.size(); ++j) {
                        def ampC = amp[j]
                        ampC <<= Conjugate & Reverse[Matrix1, Matrix2] & stp.conjugateSpinors
                        ampC <<= { expr -> (expr.indices.free.si % expr.indices.free.si.inverted) >> expr } as Transformation

                        def amp2 = amp[i] * ampC
                        def indexless = amp2.indexlessSubProduct
                        def tensor = amp2.dataSubProduct

                        tensor <<= stp.epsSum & stp.uTrace & mndlst & dTraceSimplify &
                                fullSimplify & stp.massesSubs & stp.uSimplify
                        sum << indexless * tensor
                    }

                def amp2 = sum.build()
                amp2 <<= stp.massesSubs
                amp2 <<= 'u = 2*mc**2 - s - t'.t

                assert stp.wolframFactorTr >> amp2 == 0.t
            }
        }
    }

    /**
     * In case of axial meson, check Landau-Yang theorem
     */
    @Test
    public void testEffectiveVertexAxial() {
        use(Redberry) {
            def stp = new Setup(true)
            def vertices = stp.effectiveQuarkoniaVertices()
            def vertex = vertices['axial']
            assert ({
                def temp = vertex >> 'eps1^a[h1]*eps2^b[h2]*Aaxial_{aA bB}[bottom, k1_i, k2_i]'.t
                println TensorUtils.getAllDiffSimpleTensors(temp)
                ('p_a[bottom] = k1_a + k2_a'.t
                        & stp.fullSimplify
                        & 'm[bottom] = mb'.t.hold & 's = (2*mb)**2'.t
                        & Factor
                ) >> temp
            }()) == 0.t
        }
    }

    @Test
    public void testPolarizations1() throws Exception {
        //gluon polarizations
        use(Redberry) {
            def stp = new Setup(true);
            for (def g in [1, 2]) {
                def sum = "eps${g}_a[1] * eps${g}_b[-1] + eps${g}_b[1] * eps${g}_a[-1]".t
                sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s -t'.t
                assert '2*s**(-1)*k2_{a}*k1_{b}+2*s**(-1)*k2_{b}*k1_{a}-g_{ab}'.t == (stp.mFactor >> sum)
            }
        }
    }


    @Test
    public void testPolarizations1cc() throws Exception {
        //gluon polarizations
        use(Redberry) {
            def stp = new Setup(false, true);
            for (def g in [1, 2]) {
                def sum = "eps${g}_a[1] * eps${g}_b[-1] + eps${g}_b[1] * eps${g}_a[-1]".t
                sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & stp.mFactor
                assert '2*s**(-1)*k2_{a}*k1_{b}+2*s**(-1)*k2_{b}*k1_{a}-g_{ab}'.t == (stp.mFactor >> sum)
            }
        }
    }

    @Test
    public void testPolarizations2() throws Exception {
        //axial mesons polarizations
        use(Redberry) {
            def stp = new Setup(true);
            for (def fl in ['charm', 'bottom']) {
                def sum = "eps_a[$fl, 1] * eps_b[$fl, -1] + eps_a[$fl, 0] * eps_b[$fl, 0] + eps_b[$fl, 1] * eps_a[$fl, -1]".t
                sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.mFactor

                def expected = "eps_a[h[$fl]] * eps_b[h[$fl]]".t
                expected <<= stp.epsSum & stp.massesSubs
                assert (sum - expected) == 0.t
            }
        }
    }

    @Test
    public void testPolarizations2cc() throws Exception {
        //axial mesons polarizations
        use(Redberry) {
            def stp = new Setup(false, true);
            def fl = 'bottom'
            def sum = "eps_a[$fl, 1] * eps_b[$fl, -1] + eps_a[$fl, 0] * eps_b[$fl, 0] + eps_b[$fl, 1] * eps_a[$fl, -1]".t
            sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & stp.mFactor

            def expected = "eps_a[h[$fl]] * eps_b[h[$fl]]".t
            expected <<= stp.epsSum & stp.massesSubs
            assert (sum - expected) == 0.t
        }
    }


    @Test
    public void testPolarizations3() throws Exception {
        //tensor mesons polarizations
        use(Redberry) {
            def stp = new Setup(true, true);
            for (def fl in ['bottom']) {
                def sum = new SumBuilder()
                def pairs = [["eps_ab[$fl, 2]", "eps_cd[$fl, -2]"].t,
                             ["eps_ab[$fl, 1]", "eps_cd[$fl, -1]"].t,
                             ["eps_ab[$fl, 0]", "eps_cd[$fl, 0]"].t,
                             ["eps_ab[$fl, -1]", "eps_cd[$fl, 1]"].t,
                             ["eps_ab[$fl, -2]", "eps_cd[$fl, 2]"].t]

                for (def pair in pairs) {
                    def a = pair[0], b = pair[1]
                    a <<= stp.polarisations; b <<= stp.polarisations;
                    for (int i = 0; i < a.size(); ++i)
                        for (int j = 0; j < b.size(); ++j) {
                            Product p = a[i] * b[j]
                            def ind = p.dataSubProduct
                            ind <<= stp.fullSimplify & stp.massesSubs

                            if (ind.class == Sum)
                                sum << FastTensors.multiplySumElementsOnFactor(ind, p.indexlessSubProduct)
                            else
                                sum << ind * p.indexlessSubProduct
                        }
                }

                //eps_ab * eps_cd
                sum = sum.build()
                println sum.size()

                sum <<= 'u = 4*mc**2 + 4*mb**2 - s - t'.t & stp.mFactor & stp.mFactor
                sum.each {
                    println it
                }
            }
        }
    }

    @Test
    public void testPolarizations3сс() throws Exception {
        //tensor mesons polarizations
        use(Redberry) {
            def stp = new Setup(false, true);
            for (def fl in ['bottom']) {
                def sum = new SumBuilder()
                def pairs = [["eps_ab[$fl, 2]", "eps_cd[$fl, -2]"].t,
                             ["eps_ab[$fl, 1]", "eps_cd[$fl, -1]"].t,
                             ["eps_ab[$fl, 0]", "eps_cd[$fl, 0]"].t,
                             ["eps_ab[$fl, -1]", "eps_cd[$fl, 1]"].t,
                             ["eps_ab[$fl, -2]", "eps_cd[$fl, 2]"].t]

                for (def pair in pairs) {
                    def a = pair[0], b = pair[1]
                    a <<= stp.polarisations; b <<= stp.polarisations;
                    for (int i = 0; i < a.size(); ++i)
                        for (int j = 0; j < b.size(); ++j) {
                            Product p = a[i] * b[j]
                            def ind = p.dataSubProduct
                            ind <<= stp.fullSimplify & stp.massesSubs

                            if (ind.class == Sum)
                                sum << FastTensors.multiplySumElementsOnFactor(ind, p.indexlessSubProduct)
                            else
                                sum << ind * p.indexlessSubProduct
                        }
                }

                //eps_ab * eps_cd
                sum = sum.build()
                println sum.size()

                sum <<= stp.mFactor & stp.mFactor
                sum.each {
                    println it
                }
            }
        }
    }

    @Test
    public void testPolarizationsTensorSum() throws Exception {
        //tensor mesons polarizations
        use(Redberry) {
            def stp = new Setup(true, true);
            for (def fl in ['charm', 'bottom']) {
                def pairs = [["eps_ab[$fl, 2]", "eps^ab[$fl, -2]"].t,
                             ["eps_ab[$fl, 1]", "eps^ab[$fl, -1]"].t,
                             ["eps_ab[$fl, 0]", "eps^ab[$fl, 0]"].t,
                             ["eps_ab[$fl, -1]", "eps^ab[$fl, 1]"].t,
                             ["eps_ab[$fl, -2]", "eps^ab[$fl, 2]"].t]

                for (def pair in pairs) {
                    def sum = new SumBuilder()
                    def a = pair[0], b = pair[1]
                    a <<= stp.polarisations; b <<= stp.polarisations;
                    for (int i = 0; i < a.size(); ++i)
                        for (int j = 0; j < b.size(); ++j) {
                            Product p = a[i] * b[j]
                            def ind = p.dataSubProduct
                            ind <<= stp.fullSimplify & stp.massesSubs
                            sum << p.indexlessSubProduct * ind
                        }

                    def res = sum.build()
                    res <<= 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.wolframFactorTr
                    assert res == 1.t
                }
            }
        }
    }


    @Test
    public void testPolarizationsTensorSum2() throws Exception {
        //tensor mesons polarizations
        use(Redberry) {
            def stp = new Setup(true, true);
            for (def fl in ['charm', 'bottom']) {

                for (def l1 in [-2, -1, 0, 1, 2])
                    for (def l2 in [-2, -1, 0, 1, 2]) {
                        def sum = new SumBuilder()
                        def a = "eps_ab[$fl, $l1]".t, b = "eps^ab[$fl, $l2]".t
                        a <<= stp.polarisations; b <<= stp.polarisations;
                        for (int i = 0; i < a.size(); ++i)
                            for (int j = 0; j < b.size(); ++j) {
                                Product p = a[i] * b[j]
                                def ind = p.dataSubProduct
                                ind <<= stp.fullSimplify & stp.massesSubs
                                sum << p.indexlessSubProduct * ind
                            }
                        def res = sum.build()
                        res <<= 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.wolframFactorTr
                        if (l1 == -l2)
                            assert res == 1.t
                        else
                            assert res == 0.t
                    }
            }
        }
    }

    @Test
    public void testPolarizationsTensorSum2cc() throws Exception {
        //tensor mesons polarizations
        use(Redberry) {
            def stp = new Setup(false, true);
            def fl = 'bottom'

            for (def l1 in [-1, 0, 1, 2, -2])
                for (def l2 in [-1, 0, 1, 2, -2]) {
                    println "$l1 $l2"
                    def sum = new SumBuilder()
                    def a = "eps_ab[$fl, $l1]".t, b = "eps^ab[$fl, $l2]".t
                    a <<= stp.polarisations; b <<= stp.polarisations;
                    for (int i = 0; i < a.size(); ++i)
                        for (int j = 0; j < b.size(); ++j) {
                            Product p = a[i] * b[j]
                            def ind = p.dataSubProduct
                            ind <<= stp.fullSimplify & stp.massesSubs
                            sum << p.indexlessSubProduct * ind
                        }
                    def res = sum.build()
                    res <<= stp.wolframFactorTr
                    if (l1 == -l2)
                        assert res == 1.t
                    else
                        assert res == 0.t
                }
        }
    }

    @Test
    public void testPolarizationsTensorSym() throws Exception {
        //tensor mesons polarizations
        use(Redberry) {
            def stp = new Setup(true, true);
            for (def fl in ['charm', 'bottom']) {

                for (def l in [-2, -1, 0, 1, 2]) {
                    def a = "eps_ab[$fl, $l] - eps_ba[$fl, $l]".t
                    a <<= stp.polarisations & stp.mFactor
                    assert a == 0.t

                    a = "eps_ab[$fl, $l] * p^a[$fl]".t
                    a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs & stp.mFactor
                    assert a == 0.t

                    a = "eps_a^a[$fl, $l]".t
                    a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs
                    a <<= 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.wolframFactorTr
                    assert a == 0.t
                }
            }
        }
    }

    @Test
    public void testPolarizationsTensorSymcc() throws Exception {
        //tensor mesons polarizations
        use(Redberry) {
            def stp = new Setup(true, true);
            def fl = 'bottom'
            for (def l in [-2, -1, 0, 1, 2]) {
                def a = "eps_ab[$fl, $l] - eps_ba[$fl, $l]".t
                a <<= stp.polarisations & stp.mFactor
                assert a == 0.t

                a = "eps_ab[$fl, $l] * p^a[$fl]".t
                a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs & stp.mFactor
                assert a == 0.t

                a = "eps_a^a[$fl, $l]".t
                a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs
                a <<= stp.wolframFactorTr
                assert a == 0.t
            }
        }
    }
}
